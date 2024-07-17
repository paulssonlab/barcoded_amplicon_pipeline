from collections import defaultdict
from functools import partial

import polars as pl
import pyarrow as pa
from cytoolz import compose
from gfapy import Gfa

from paulssonlab.sequencing.align import pairwise_align
from paulssonlab.sequencing.cigar import (
    OP_TO_COLUMN_NAME,
    _parse_variant,
    cut_cigar,
    decode_cigar,
)
from paulssonlab.sequencing.gfa import assemble_seq_from_path, gfa_name_mapping
from paulssonlab.sequencing.io import iter_bam_and_gaf


def categorical_list_hash(col):
    return col.cast(pl.List(pl.String)).list.eval(pl.element().hash()).hash()


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
):
    (
        path_filter,
        reverse_segments,
        reverse_path_mapping,
    ) = _segments_for_normalize_path(forward_segments)
    df = (
        df.rename({"path": "full_path"})
        .with_columns(
            pl.col("full_path")
            .list.set_intersection(pl.lit(path_filter, dtype=pl.List(pl.Categorical)))
            .alias("_path"),
        )
        .with_columns(
            pl.col("_path")
            .list.reverse()
            .list.eval(
                pl.element().replace(
                    reverse_path_mapping, default=None, return_dtype=pl.Categorical
                )
            )
            .alias("_path_reversed"),
            (
                pl.col("full_path")
                .list.set_intersection(
                    pl.lit(forward_segments, dtype=pl.List(pl.Categorical))
                )
                .list.len()
                > 0
            ).alias("_is_forward"),
            (
                pl.col("full_path")
                .list.set_intersection(
                    pl.lit(reverse_segments, dtype=pl.List(pl.Categorical))
                )
                .list.len()
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
        .select(
            pl.all().exclude("_path", "_path_reversed", "_is_forward", "_is_reverse")
        )
    )
    return df


def flag_end_to_end(df, endpoints, path_column="path"):
    if len(endpoints) != 2:
        raise ValueError("endpoints should contain two lists of segment names")
    if not endpoints[0] or not endpoints[1]:
        raise ValueError(f"both lists of endpoints must be non-empty: {endpoints}")
    return df.with_columns(
        end_to_end=(
            pl.col(path_column)
            .list.set_intersection(pl.lit(endpoints[0], dtype=pl.List(pl.Categorical)))
            .list.len()
            > 0
        )
        & (
            pl.col(path_column)
            .list.set_intersection(pl.lit(endpoints[1], dtype=pl.List(pl.Categorical)))
            .list.len()
            > 0
        )
    )


def flag_valid_ont_duplex_reads(df):
    # output makes the most sense if input has only
    # end-to-end (barcode) alignments, nonduplicated alignments
    df_input = df.with_columns(
        pl.col("name").str.split(";").alias("_parent_names"),
    )
    df_with_parents = (
        df_input.filter(pl.col("dx") == 1)
        .select(
            pl.col("name"),
            pl.col("path").alias("_duplex_path"),
            pl.col("_parent_names").alias("_parent_name"),
        )
        .explode("_parent_name")
        .join(
            df.select(
                pl.col("name"),
                pl.col("reverse_complement"),
                pl.col("path").alias("_parent_path"),
            ),
            how="left",
            left_on=pl.col("_parent_name"),
            right_on=pl.col("name"),
        )
    )
    df_duplex_is_valid = (
        df_with_parents.with_columns(
            (pl.col("_duplex_path") == pl.col("_parent_path")).alias("_path_matches")
        )
        .group_by("name")
        .agg(
            pl.col("_path_matches").all(),
            pl.len().alias("_group_size"),
            pl.col("reverse_complement").sum().alias("_num_rc"),
            pl.col("_parent_name"),
        )
        .with_columns(
            pl.col("_path_matches")
            .and_(pl.col("_group_size") == 2, pl.col("_num_rc") == 1)
            .alias("_duplex_is_valid")
        )
    )
    df_simplex_is_valid = df_duplex_is_valid.select(
        pl.col("_duplex_is_valid").not_().alias("_simplex_is_valid"),
        pl.col("_parent_name"),
    ).explode("_parent_name")
    df_valid = (
        df.join(
            df_duplex_is_valid.select(pl.col("name", "_duplex_is_valid")),
            how="left",
            on="name",
        )
        .join(df_simplex_is_valid, how="left", left_on="name", right_on="_parent_name")
        .with_columns(
            pl.when(pl.col("dx") == 1)
            .then(pl.col("_duplex_is_valid").fill_null(False))
            .otherwise(pl.col("_simplex_is_valid").fill_null(True))
            .alias("is_valid")
        )
        .select(pl.all().exclude("_duplex_is_valid", "_simplex_is_valid"))
    )
    return df_valid


def compute_divergence(segment, col_getter=pl.col):
    return pl.sum_horizontal(
        col_getter(f"{segment}|{type_}").fill_null(strategy="zero")
        for type_ in ("mismatches", "insertions", "deletions")
    ) / pl.sum_horizontal(
        col_getter(f"{segment}|{type_}").fill_null(strategy="zero")
        for type_ in ("matches", "mismatches", "insertions", "deletions")
    )


def compute_divergences(df, segments, struct_name=None):
    if struct_name is not None:
        col_getter = pl.col(struct_name).struct.field
    else:
        col_getter = pl.col
    cols = {
        f"{segment}|divergence": compute_divergence(
            segment,
            col_getter=col_getter,
        )
        for segment in segments
    }
    if struct_name is not None:
        df = df.with_columns(pl.col(struct_name).struct.with_fields(**cols))
        divergence_cols = [
            pl.col(struct_name).struct.field(f"{segment}|divergence")
            for segment in segments
        ]
    else:
        df = df.with_columns(**cols)
        divergence_cols = [pl.col(f"{segment}|divergence") for segment in segments]
    df = df.with_columns(max_divergence=pl.max_horizontal(*divergence_cols))
    return df


def unique_segments(df, col):
    if isinstance(col, str):
        col = pl.col(col)
    segments = df.select(
        col.explode()
        .cast(pl.String)
        .str.slice(1, None)
        .str.split("=")
        .list.get(0)
        .unique()
        .alias("segment")
    )
    if hasattr(segments, "collect"):
        segments = segments.collect()
    return segments.to_series().to_list()


def prepare_reads(df, forward_segments, endpoints, name_to_seq, max_divergence=None):
    df = normalize_paths(df, forward_segments)
    df = df.with_columns(path_hash=categorical_list_hash(pl.col("path")))
    df = flag_end_to_end(df, endpoints)
    # TODO: is there any reason to make these column names configurable?
    df = cut_cigar_df(
        df,
        name_to_seq,
        path_column="full_path",
        cigar_column="cg",
        # sequence_column="read_seq",
        # phred_column="read_phred",
        query_start_column="query_start",
        query_end_column="query_end",
        query_length_column="query_length",
        path_start_column="path_start",
        path_end_column="path_end",
        keep_full=True,
        struct_name="extract_segments",
        cut_cigar_kwargs=dict(
            key_sep="|",
            return_indices=False,
            return_counts=True,
            return_cigars=False,
            return_variants=True,
            separate_ends=True,
        ),
    )
    # --max-divergence must be specified to prepare_reads.py for that filtering to be used
    # for filtering duplex reads when use_dorado_duplex_pairing=True
    # otherwise --max-divergence need only be passed to consensus.py when computing consensus seqs
    candidate = pl.col("is_primary_alignment") & pl.col("end_to_end")
    if max_divergence is not None:
        df = compute_divergences(
            df, [s[1:] for s in forward_segments], struct_name="extract_segments"
        )
        candidate = candidate & (pl.col("max_divergence") <= max_divergence)
    df = df.with_columns(
        is_primary_alignment=pl.col("name").is_first_distinct(),
    ).with_columns(_candidate=candidate)
    # only ONT duplex reads will have dx SAM tag
    # TODO
    # if "dx" in df.collect_schema().names():
    if "dx" in df.columns:
        df_valid_reads = flag_valid_ont_duplex_reads(df.filter(pl.col("_candidate")))
        df = pl.concat(
            [df_valid_reads, df.filter(~pl.col("_candidate"))], how="diagonal"
        )
    # TODO
    if "qs" not in df.columns and "read_phred" in df.columns:
        df = df.with_columns(qs=pl.col("read_phred").list.mean().cast(pl.Float32))
    df = df.select(pl.all().exclude("_candidate"))
    return df


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


def find_duplex_pairs(
    df,
    timedelta,
    forward_segments,
    endpoints=None,
    max_divergence=None,
    name_to_seq=None,
):
    (
        path_filter,
        _,
        reverse_path_mapping,
    ) = _segments_for_normalize_path(forward_segments)
    df = (
        df.rename({"path": "full_path"})
        .filter(pl.col("name").is_first_distinct())  # use primary alignment
        .with_columns(
            pl.col("full_path")
            .list.set_intersection(pl.lit(path_filter, dtype=pl.List(pl.Categorical)))
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
        df = flag_end_to_end(df, endpoints, path_column="_path").filter(
            pl.col("end_to_end")
        )
    if max_divergence is not None:
        df = cut_cigar_df(
            df,
            name_to_seq,
            path_column="full_path",
            cigar_column="cg",
            # sequence_column="read_seq",
            # phred_column="read_phred",
            query_start_column="query_start",
            query_end_column="query_end",
            query_length_column="query_length",
            path_start_column="path_start",
            path_end_column="path_end",
            keep_full=True,
            struct_name="extract_segments",
            cut_cigar_kwargs=dict(
                key_sep="|",
                return_indices=False,
                return_counts=True,
                return_cigars=False,
                return_variants=True,
                separate_ends=True,
            ),
        )
        df = compute_divergences(
            df, [s[1:] for s in forward_segments], struct_name="extract_segments"
        )
        df = df.filter(pl.col("max_divergence") <= max_divergence)
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
        tolerance=timedelta,
    ).filter(
        (pl.col("_path") == pl.col("_path_reversed_right"))
        & pl.col("name_right").is_not_null()
    )
    return df_joined


def compute_depth(df, over=None, prefix=None, suffix=None):
    if prefix is None:
        prefix = ""
    if suffix is None:
        suffix = ""

    def format_column(name):
        return f"{prefix}{name}{suffix}"

    if over is not None:

        def wrap_over(expr):
            return expr.over(over)

    else:

        def wrap_over(expr):
            return expr

    exprs = {format_column("depth"): wrap_over(pl.len()).alias(format_column("depth"))}
    # TODO
    # if "dx" in df.collect_schema().names():
    if "dx" in df.columns:
        exprs[format_column("duplex_depth")] = wrap_over(
            (pl.col("dx") == 1).sum().alias(format_column("duplex_depth"))
        )
    return exprs


def map_read_groups(
    df,
    func,
    *agg_exprs,
    max_group_size=None,
    input_columns=["name", "read_seq", "read_phred", "reverse_complement"],
):
    # TODO: this will be slightly faster/cleaner when polars implements top_k for structs
    # SEE: https://github.com/pola-rs/polars/issues/7043
    # and https://stackoverflow.com/questions/76596952/how-do-i-select-the-top-k-rows-of-a-python-polars-dataframe-for-each-group
    # and https://github.com/pola-rs/polars/issues/10054
    def _limit_group(expr):
        sort_columns = []
        descending = []
        # TODO
        # df_columns = df.collect_schema().names()
        df_columns = df.columns
        # PacBio: prefer high-accuracy CCS reads
        if "rq" in df_columns:
            sort_columns = ["rq", *sort_columns]
            descending = [True, *descending]
        # ONT: prefer high-accuracy reads
        # PacBio CCS uses the "qs" tag for something else, so ignore if "rq" is present
        if "qs" in df_columns and "rq" not in df_columns:
            sort_columns = ["qs", *sort_columns]
            descending = [True, *descending]
        if sort_columns:
            expr = expr.sort_by(
                sort_columns,
                descending=descending,
            )
        if max_group_size is not None:
            expr = expr.head(max_group_size)
        # we don't sort by length here because:
        # 1) consensus.poa sorts by length already
        # 2) there's a polars bug, this doesn't work for some reason
        # (maybe it doesn't like multiple sort_by's in a group_by?)
        # expr = expr.sort_by(pl.col("read_seq").str.len_bytes(), descending=True)
        return expr

    return (
        df.group_by("path")
        .agg(_limit_group(pl.all()))
        .explode(pl.all().exclude("path"))
        .group_by("path")
        .agg(
            pl.map_groups(
                pl.struct(pl.col(input_columns)),
                compose(func, lambda series: series[0].struct.unnest()),
                returns_scalar=True,
            ).alias("_func_output"),
            *agg_exprs,
        )
    ).unnest("_func_output")


def _pairwise_align_row(res, score_column, cigar_column):
    return {score_column: res["score"], cigar_column: res["cigar"]}


def _pairwise_align_rows(
    rows,
    name_to_seq=None,
    score_column=None,
    cigar_column=None,
    dtype=None,
    align_kwargs={},
):
    align_kwargs = {"cigar_as_string": True, **align_kwargs}
    paths = rows.struct.field("path")
    seqs = rows.struct.field("seq")
    return pl.Series(
        [
            _pairwise_align_row(
                pairwise_align(
                    seqs[idx],
                    assemble_seq_from_path(name_to_seq, paths[idx]),
                    **align_kwargs,
                ),
                score_column,
                cigar_column,
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
    dtype = pl.Struct({score_column: pl.Int32, cigar_column: pl.String})
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


def _slice_if_not_none(obj, idx):
    if obj is not None and idx is not None:
        return obj[idx]
    else:
        return None


def _get_field(rows, name):
    if name in rows.struct.fields:
        return rows.struct.field(name)
    else:
        return None


def _cut_cigar_rows(rows, name_to_seq=None, cut_cigar_kwargs={}, dtype=None):
    paths = rows.struct.field("path")
    cigars = rows.struct.field("cigar")
    seqs = _get_field(rows, "seq")
    phreds = _get_field(rows, "phred")
    query_start = _get_field(rows, "query_start")
    query_end = _get_field(rows, "query_end")
    query_length = _get_field(rows, "query_length")
    path_start = _get_field(rows, "path_start")
    path_end = _get_field(rows, "path_end")
    return pl.Series(
        [
            (
                cut_cigar(
                    decode_cigar(cigars[idx]),
                    paths[idx].to_list(),
                    name_to_seq,
                    sequence=_slice_if_not_none(seqs, idx),
                    phred=phreds[idx].to_numpy() if phreds is not None else None,
                    query_start=_slice_if_not_none(query_start, idx),
                    query_end=_slice_if_not_none(query_end, idx),
                    query_length=_slice_if_not_none(query_length, idx),
                    path_start=_slice_if_not_none(path_start, idx),
                    path_end=_slice_if_not_none(path_end, idx),
                    cigar_as_string=True,
                    **cut_cigar_kwargs,
                )
                if _slice_if_not_none(cigars, idx) is not None
                and _slice_if_not_none(paths, idx) is not None
                else None
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
            segments_to_variant_type[k] = pl.String
    for segment_name in segment_names:
        if segment_name in full_to_segment_name:
            segment_name = full_to_segment_name[segment_name]
            if return_variants and variants[segment_name]:
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
            dtype[f"{segment_name}{key_sep}cigar"] = pl.String
    return pl.Struct(dtype)


def _include_column(columns, exclude_columns, source_columns, source_name, dest_name):
    if source_name is not None and source_name in source_columns:
        columns[dest_name] = source_name
        if exclude_columns is not None:
            exclude_columns.append(source_name)


def cut_cigar_df(
    df,
    name_to_seq,
    path_column=None,
    cigar_column=None,
    sequence_column=None,
    phred_column=None,
    query_start_column=None,
    query_end_column=None,
    query_length_column=None,
    path_start_column=None,
    path_end_column=None,
    struct_name=None,
    keep_full=False,
    cut_cigar_kwargs={},
):
    if isinstance(name_to_seq, Gfa):
        name_to_seq = gfa_name_mapping(name_to_seq)
    exclude_columns = []
    # TODO
    # df_columns = df.collect_schema().names()
    df_columns = df.columns
    if path_column not in df_columns:
        raise ValueError(f"missing column {path_column}")
    if cigar_column not in df_columns:
        raise ValueError(f"missing column {cigar_column}")
    else:
        exclude_columns.append(cigar_column)
    struct = dict(path=path_column, cigar=cigar_column)
    _include_column(struct, exclude_columns, df_columns, sequence_column, "seq")
    _include_column(struct, exclude_columns, df_columns, phred_column, "phred")
    _include_column(struct, None, df_columns, query_start_column, "query_start")
    _include_column(struct, None, df_columns, query_end_column, "query_end")
    _include_column(struct, None, df_columns, query_length_column, "query_length")
    _include_column(struct, None, df_columns, path_start_column, "path_start")
    _include_column(struct, None, df_columns, path_end_column, "path_end")
    if keep_full:
        exclude_columns = []
    # TODO
    # df_schema = df.collect_schema()
    df_schema = df.schema
    dtype = _cut_cigar_dtype(
        [s[1:] for s in name_to_seq.keys()],
        sequence_dtype=df_schema.get(sequence_column),
        phred_dtype=df_schema.get(phred_column),
        **cut_cigar_kwargs,
    )
    df = df.select(
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
        .alias("_func_output" if struct_name is None else struct_name),
    )
    if struct_name is None:
        df = df.unnest("_func_output")
    return df
