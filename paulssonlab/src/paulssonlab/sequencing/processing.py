from functools import partial

import polars as pl

from paulssonlab.sequencing.align import pairwise_align
from paulssonlab.sequencing.cigar import cut_cigar, decode_cigar
from paulssonlab.sequencing.gfa import assemble_seq_from_path, gfa_name_mapping


# inspired by https://github.com/pola-rs/polars/issues/8205#issuecomment-1508532138
def join_dfs(
    left,
    right,
    on=None,
    left_prefix=None,
    left_suffix=None,
    right_prefix=None,
    right_suffix=None,
    **kwargs,
):
    if on is None:
        raise ValueError("must specify on")
    if isinstance(on, str):
        on = [on]
    if left_prefix is None:
        left_prefix = ""
    if left_suffix is None:
        left_suffix = ""
    if right_prefix is None:
        right_prefix = ""
    if right_suffix is None:
        right_suffix = ""
    common_columns = set(left.columns) & set(right.columns) - set(on)
    left_renamed = left.rename(
        {
            k: f"{left_prefix}{k}{left_suffix}" if k in common_columns else k
            for k in left.columns
        }
    )
    right_renamed = right.rename(
        {
            k: f"{right_prefix}{k}{right_suffix}" if k in common_columns else k
            for k in right.columns
        }
    )
    return left_renamed.join(right_renamed, on=on, **kwargs)


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


def normalize_path(
    df,
    forward_segments,
    keep_unaligned=False,
    endpoints=None,
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


def map_read_groups(df, func):
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
                pl.struct(
                    pl.col("name", "read_seq", "read_phred", "reverse_complement")
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


def cut_cigar_exc(*args, **kwargs):
    try:
        return cut_cigar(*args, **kwargs)
    except:
        print("ARGS", args)
        print("KWARGS", kwargs)
        print()
        print()
        raise


def _cut_cigar_rows(rows, name_to_seq=None, cut_cigar_kwargs={}):
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
            cut_cigar_exc(
                decode_cigar(cigars[idx]),
                paths[idx].to_list(),
                name_to_seq,
                sequence=seqs[idx] if seqs is not None else None,
                phred=phreds[idx].to_numpy() if phreds is not None else None,
                cigar_as_string=True,
                **cut_cigar_kwargs,
            )
            for idx in range(len(rows))
        ]
    )


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
    if path_column not in df.columns:
        raise ValueError(f"missing column {path_column}")
    if cigar_column not in df.columns:
        raise ValueError(f"missing column {cigar_column}")
    exclude_columns = []
    struct = dict(path=path_column, cigar=cigar_column)
    if sequence_column and sequence_column in df.columns:
        exclude_columns.append(sequence_column)
        struct["seq"] = sequence_column
    if phred_column and phred_column in df.columns:
        struct["phred"] = phred_column
        exclude_columns.append(phred_column)
    exclude_columns = []
    if keep_full:
        exclude_columns = []
    return (
        df.head(1000)
        .select(
            pl.all().exclude(exclude_columns),
            pl.struct(**struct)
            # TODO: dtype depends on cut_cigar_kwargs in complicated ways
            # so we don't supply it
            # polars docs says this is an error but still supported (??)
            # in any case, it works for now
            .map_batches(
                partial(
                    _cut_cigar_rows,
                    name_to_seq=name_to_seq,
                    cut_cigar_kwargs=cut_cigar_kwargs,
                ),
            ).alias("_func_output"),
        )
        .unnest("_func_output")
    )
