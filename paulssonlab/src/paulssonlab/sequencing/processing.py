from collections import defaultdict
from functools import partial

import polars as pl
import pyarrow as pa

from paulssonlab.sequencing.align import pairwise_align
from paulssonlab.sequencing.cigar import (
    OP_TO_COLUMN_NAME,
    _parse_variant,
    cut_cigar,
    decode_cigar,
)
from paulssonlab.sequencing.gfa import assemble_seq_from_path, gfa_name_mapping
from paulssonlab.sequencing.io import iter_bam_and_gaf


def _reverse_segment(s):
    if s[0] == ">":
        return f"<{s[1:]}"
    elif s[0] == "<":
        return f">{s[1:]}"
    else:
        raise ValueError(f"segment name should start with < or >: {s}")


def _segments_for_normalize_path(forward_segments):
    segment_names = [s[1:] for s in forward_segments]
    path_filter = [f"{o}{s}" for s in segment_names for o in "><"]
    reverse_segments = [_reverse_segment(s) for s in forward_segments for o in "><"]
    reverse_path_mapping = {
        f"{o[0]}{s}": f"{o[1]}{s}" for s in segment_names for o in ["<>", "><"]
    }
    return path_filter, reverse_segments, reverse_path_mapping


def normalize_paths(
    df,
    forward_segments,
    keep_unaligned=False,
    endpoints=None,
    hash_paths=True,
):
    if endpoints:
        if len(endpoints) != 2:
            raise ValueError("endpoints should contain two lists of segment names")
        if not endpoints[0] or not endpoints[1]:
            raise ValueError(f"both lists of endpoints must be non-empty: {endpoints}")
    (
        path_filter,
        reverse_segments,
        reverse_path_mapping,
    ) = _segments_for_normalize_path(forward_segments)
    df = (
        df.rename({"path": "full_path"})
        # .is_first_distinct() picks better (primary) alignment instead of alternate alignment
        # if input is in GAF order
        # .filter(pl.col("name").is_first_distinct())
        # TODO: polars bug, is_first_distinct() returns empty dataframe
        # when used in a more complicated lazy query,
        # exact cause unknown
        .filter(pl.col("name").is_duplicated().not_())
        .with_columns(
            pl.col("full_path")
            .list.set_intersection(path_filter)
            # TODO: waiting on https://github.com/pola-rs/polars/issues/11735
            # to keep path columns as categorical
            .cast(pl.List(pl.Categorical))
            .alias("_path"),
        )
        .with_columns(
            pl.col("_path")
            .list.reverse()
            .list.eval(pl.element().map_dict(reverse_path_mapping))
            .alias("_path_reversed"),
            (
                pl.col("full_path").list.set_intersection(forward_segments).list.len()
                > 0
            ).alias("_is_forward"),
            (
                pl.col("full_path").list.set_intersection(reverse_segments).list.len()
                > 0
            ).alias("_is_reverse"),
        )
        .with_columns(
            pl.when(pl.col("_is_forward") & pl.col("_is_reverse").not_())
            .then(False)
            .when(pl.col("_is_forward").not_() & pl.col("_is_reverse"))
            .then(True)
            .otherwise(None)
            .alias("reverse_complement"),
        )
        .with_columns(
            pl.when(pl.col("reverse_complement") == False)
            .then(pl.col("_path"))
            .when(pl.col("reverse_complement") == True)
            .then(pl.col("_path_reversed"))
            .alias("path")
        )
    )
    if hash_paths:
        df = df.with_columns(pl.col("path").hash().alias("path_hash"))
    if not keep_unaligned:
        df = df.filter(pl.col("path").is_not_null())
    if endpoints:
        df = df.filter(
            (pl.col("path").list.set_intersection(endpoints[0]).list.len() > 0)
            & (pl.col("path").list.set_intersection(endpoints[1]).list.len() > 0)
        )
    return df


def identify_usable_reads(df):
    df_input = df.with_columns(
        pl.col("name").str.split(";").alias("parent_names")
    ).with_columns(
        pl.when(pl.col("parent_names").list.len() == 2)
        .then(True)
        .otherwise(False)
        .alias("is_duplex")
    )
    df_with_parents = (
        df_input.filter(pl.col("is_duplex"))
        .select(pl.col("name"), pl.col("parent_names").alias("_parent_name"))
        .explode("_parent_name")
        .join(
            df.select(pl.col("name", "path")),
            how="left",
            left_on=pl.col("_parent_name"),
            right_on=pl.col("name"),
        )
    )
    df_duplex_paths_match = (
        df_with_parents.with_columns(
            pl.col("path").first().over("name").alias("_path_first")
        )
        .with_columns(
            (pl.col("path") == pl.col("_path_first")).alias("_path_matches_first")
        )
        .group_by("name")
        .agg(
            pl.col("_path_matches_first").all().alias("_duplex_paths_match"),
            pl.col("_parent_name"),
        )
    )
    df_simplex_paths_match = df_duplex_paths_match.select(
        pl.col("_duplex_paths_match").alias("_child_duplex_paths_match"),
        pl.col("_parent_name"),
    ).explode("_parent_name")
    df_usable = (
        df_input.join(
            df_duplex_paths_match.select(pl.col("name", "_duplex_paths_match")),
            how="left",
            on="name",
        )
        .join(
            df_simplex_paths_match, how="left", left_on="name", right_on="_parent_name"
        )
        .with_columns(
            pl.when(pl.col("is_duplex"))
            .then(pl.col("_duplex_paths_match").fill_null(False))
            .otherwise(pl.col("_child_duplex_paths_match").not_().fill_null(True))
            .alias("usable")
        )
    )
    return df_usable


def compute_depth(df):
    return df.with_columns(
        pl.col("is_duplex").sum().over("path").alias("duplex_depth"),
        pl.count().over("path").alias("depth"),
    ).with_columns(
        (pl.col("depth") - pl.col("duplex_depth")).alias("simplex_depth"),
    )


def load_pairing_data(bam_filename, gaf_filename):
    batches = []
    for batch in iter_bam_and_gaf(bam_filename, gaf_filename, include_unaligned=False):
        # TODO: could do this processing more cleanly in polars
        # but wasn't sure how to do pl.LazyFrame computations on an iterator of RecordBatches
        batch = batch.select(["name", "path", "qs", "mx", "ch", "st", "du", "pi"])
        schema = batch.schema
        schema = schema.set(
            schema.names.index("st"), pa.field("st", pa.timestamp("ms", "UTC"))
        )
        batch = pa.RecordBatch.from_pydict(
            {
                **batch.to_pydict(),
                "st": batch.column("st").cast(pa.timestamp("ms", "UTC")),
            },
            schema=schema,
        )
        batches.append(batch)
    return pl.from_arrow(batches)


def find_duplex_pairs(df, forward_segments, endpoints=None):
    (
        path_filter,
        _,
        reverse_path_mapping,
    ) = _segments_for_normalize_path(forward_segments)
    df = (
        df.rename({"path": "full_path"})
        # .is_first_distinct() picks better (primary) alignment instead of alternate alignment
        # if input is in GAF order
        # .filter(pl.col("name").is_first_distinct())
        # TODO: polars bug, is_first_distinct() returns empty dataframe
        # when used in a more complicated lazy query,
        # exact cause unknown
        .filter(pl.col("name").is_first_distinct())
        .with_columns(
            pl.col("full_path")
            .list.set_intersection(path_filter)
            # TODO: waiting on https://github.com/pola-rs/polars/issues/11735
            # to keep path columns as categorical
            .cast(pl.List(pl.Categorical))
            .alias("_path"),
        )
        .with_columns(
            pl.col("_path")
            .list.reverse()
            .list.eval(pl.element().replace(reverse_path_mapping, default=None))
            .alias("_path_reversed"),
        )
    )
    if endpoints:
        df = df.filter(
            (pl.col("_path").list.set_intersection(endpoints[0]).list.len() > 0)
            & (pl.col("_path").list.set_intersection(endpoints[1]).list.len() > 0)
        )
    df_sorted = df.sort("st")
    df_endtime = df.with_columns(
        (pl.col("st") + pl.duration(seconds="du")).alias("_endtime")
    ).sort("_endtime")
    df_joined = df_endtime.join_asof(
        df_sorted,
        left_on="_endtime",
        right_on="st",
        by=["ch", "mx"],
        strategy="forward",
        tolerance="10s",
    ).filter(
        (pl.col("_path") == pl.col("_path_reversed_right"))
        & pl.col("name_right").is_not_null()
    )
    return df_joined


def map_read_groups(df, func, max_group_size=None):
    if max_group_size is None:

        def _limit_group(expr):
            return expr

    else:
        # prefer duplex reads, longer reads
        def _limit_group(expr):
            return expr.sort_by(
                ["is_duplex", pl.col("read_seq").str.len_bytes()],
                descending=[True, True],
            ).head(max_group_size)

    return (
        df.select(
            pl.col(
                "name",
                "is_duplex",
                "read_seq",
                "read_phred",
                "reverse_complement",
                "path",
                "depth",
                "simplex_depth",
                "duplex_depth",
            )
        )
        .group_by("path")
        .agg(
            pl.map_groups(
                _limit_group(
                    pl.struct(
                        pl.col("name", "read_seq", "read_phred", "reverse_complement")
                    )
                ),
                func,
                returns_scalar=True,
            ).alias("_func_output"),
            pl.col("path", "depth", "simplex_depth", "duplex_depth").first(),
        )
    ).unnest("_func_output")


def _pairwise_align_rows(
    rows,
    name_to_seq=None,
    score_column=None,
    cigar_column=None,
    dtype=None,
    align_kwargs={},
):
    paths = rows.struct.field("path")
    seqs = rows.struct.field("seq")
    return pl.Series(
        [
            dict(
                zip(
                    (score_column, cigar_column),
                    pairwise_align(
                        seqs[idx],
                        assemble_seq_from_path(name_to_seq, paths[idx]),
                        **align_kwargs,
                    ),
                )
            )
            for idx in range(len(rows))
        ],
        dtype=dtype,
    )


def pairwise_align_df_to_path(
    df,
    gfa,
    path_column=None,
    sequence_column=None,
    score_column=None,
    cigar_column=None,
    align_kwargs={},
):
    name_to_seq = gfa_name_mapping(gfa)
    dtype = pl.Struct({score_column: pl.Int32, cigar_column: pl.Utf8})
    return df.select(
        pl.all(),
        pl.struct(
            path=path_column,
            seq=sequence_column,
        )
        .map_batches(
            partial(
                _pairwise_align_rows,
                name_to_seq=name_to_seq,
                score_column=score_column,
                cigar_column=cigar_column,
                dtype=dtype,
                align_kwargs=align_kwargs,
            ),
            return_dtype=dtype,
        )
        .alias("_func_output"),
    ).unnest("_func_output")


def _cut_cigar_rows(rows, name_to_seq=None, cut_cigar_kwargs={}, dtype=None):
    paths = rows.struct.field("path")
    cigars = rows.struct.field("cigar")
    if "seq" in rows.struct.fields:
        seqs = rows.struct.field("seq")
    else:
        seqs = None
    if "phred" in rows.struct.fields:
        phreds = rows.struct.field("phred")
    else:
        phreds = None
    return pl.Series(
        [
            cut_cigar(
                decode_cigar(cigars[idx]),
                paths[idx].to_list(),
                name_to_seq,
                sequence=seqs[idx] if seqs is not None else None,
                phred=phreds[idx].to_numpy() if phreds is not None else None,
                cigar_as_string=True,
                **cut_cigar_kwargs,
            )
            for idx in range(len(rows))
        ],
        dtype=dtype,
    )


def _cut_cigar_dtype(
    segment_names,
    sequence_dtype=None,
    phred_dtype=None,
    segments=None,
    variant_sep="=",
    key_sep="|",
    return_indices=True,
    return_counts=True,
    return_cigars=True,
    cigar_as_string=True,
    return_variants=True,
    separate_ends=True,
):
    if not cigar_as_string:
        raise ValueError("cigar_as_string=False is not supported")
    if key_sep is None:
        raise ValueError("key_sep=None is not supported")
    dtype = {}
    if separate_ends:
        segment_names = ["upstream", *segment_names, "downstream"]
    variants = defaultdict(list)
    full_to_segment_name = {}
    for segment_name in segment_names:
        if variant_sep in segment_name:
            full_segment_name = segment_name
            sep_idx = segment_name.index(variant_sep)
            segment_name, variant_name = (
                segment_name[:sep_idx],
                _parse_variant(segment_name[sep_idx + 1 :]),
            )
            full_to_segment_name[full_segment_name] = segment_name
            variants[segment_name].append(variant_name)
    segments_to_variant_type = {}
    for k in variants:
        if all(isinstance(v, int) for v in variants[k]):
            segments_to_variant_type[k] = pl.Int32
        else:
            segments_to_variant_type[k] = pl.Utf8
    if return_variants:
        dtype[f"{segment_name}{key_sep}variant"] = pl.Utf8
    for segment_name in segment_names:
        if segment_name in full_to_segment_name:
            segment_name = full_to_segment_name[segment_name]
            if return_variants:
                dtype[f"{segment_name}{key_sep}variant"] = segments_to_variant_type[
                    segment_name
                ]
        if segments is not None and segment_name not in segments:
            continue
        if sequence_dtype is not None:
            dtype[f"{segment_name}{key_sep}seq"] = sequence_dtype
        if phred_dtype is not None:
            dtype[f"{segment_name}{key_sep}phred"] = phred_dtype
        if return_indices:
            dtype[f"{segment_name}{key_sep}start"] = pl.Int32
            dtype[f"{segment_name}{key_sep}end"] = pl.Int32
            dtype[f"{segment_name}{key_sep}reverse_complement"] = pl.Boolean
        if return_counts:
            for op in OP_TO_COLUMN_NAME.values():
                dtype[f"{segment_name}{key_sep}{op}"] = pl.Int32
        if return_cigars:
            dtype[f"{segment_name}{key_sep}cigar"] = pl.Utf8
    return pl.Struct(dtype)


def cut_cigar_df(
    df,
    gfa,
    path_column=None,
    cigar_column=None,
    sequence_column=None,
    phred_column=None,
    keep_full=False,
    cut_cigar_kwargs={},
):
    name_to_seq = gfa_name_mapping(gfa)
    exclude_columns = []
    if path_column not in df.columns:
        raise ValueError(f"missing column {path_column}")
    if cigar_column not in df.columns:
        raise ValueError(f"missing column {cigar_column}")
    else:
        exclude_columns.append(cigar_column)
    struct = dict(path=path_column, cigar=cigar_column)
    if sequence_column and sequence_column in df.columns:
        exclude_columns.append(sequence_column)
        struct["seq"] = sequence_column
    if phred_column and phred_column in df.columns:
        struct["phred"] = phred_column
        exclude_columns.append(phred_column)
    if keep_full:
        exclude_columns = []
    dtype = _cut_cigar_dtype(
        [s[1:] for s in name_to_seq.keys()],
        sequence_dtype=df.schema.get(sequence_column),
        phred_dtype=df.schema.get(phred_column),
        **cut_cigar_kwargs,
    )
    return df.select(
        pl.all().exclude(exclude_columns),
        pl.struct(**struct)
        .map_batches(
            partial(
                _cut_cigar_rows,
                name_to_seq=name_to_seq,
                cut_cigar_kwargs=cut_cigar_kwargs,
                dtype=dtype,
            ),
            return_dtype=dtype,
        )
        .alias("_func_output"),
    ).unnest("_func_output")
