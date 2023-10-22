import polars as pl

from paulssonlab.sequencing.gfa import gfa_endpoints


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
    reverse_segments = [[_reverse_segment(s) for s in segment_names] for o in "><"]
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
    if endpoints and len(endpoints) != 2:
        raise ValueError("endpoints should contain two iterables of segment names")
    (
        path_filter,
        reverse_segments,
        reverse_path_mapping,
    ) = _segments_for_normalize_path(forward_segments)
    df = (
        df.rename({"path": "full_path"})
        # .is_first_distinct() picks better (primary) alignment instead of alternate alignment
        # if input is in GAF order
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
            df_path.select(pl.col("name", "path")),
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


def group_by_path(df):
    return (
        df.select(
            pl.col(
                "name",
                "is_duplex",
                "read_seq",
                "read_phred",
                "reverse_complement",
                "path",
                "full_path",
            )
        )
        .group_by("path")
        .agg(
            pl.col("name", "read_seq", "read_phred", "reverse_complement"),
            pl.col("is_duplex").sum().alias("duplex_depth"),
            pl.count().alias("depth"),
        )
        .with_columns((pl.col("depth") - pl.col("duplex_depth")).alias("simplex_depth"))
    )
