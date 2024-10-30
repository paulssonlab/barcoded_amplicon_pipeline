#!/usr/bin/env python
import sys
from functools import partial
from pathlib import Path

import click
import polars as pl

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.consensus import get_consensus_group_by
from paulssonlab.sequencing.gfa import filter_gfa_options, filter_segments
from paulssonlab.sequencing.io import write_fastx
from paulssonlab.sequencing.processing import (
    categorical_list_hash,
    compute_depth,
    compute_divergences,
    map_read_groups,
    unique_segments,
)
from paulssonlab.sequencing.util import detect_format
from paulssonlab.util.cli import parse_kv, split_delimited_list

DEFAULT_FLAG_COLUMNS = ["is_primary_alignment", "end_to_end", "is_valid"]

LIMIT_DEPTH_COLUMNS = ["grouping_segments", "read_seq", "read_phred"]


def limit_group_content(df, columns, group_by, size, fill_value=None):
    # if we are not computing consensuses, we still include more than limit_depth number of reads
    # but set their read_seq and read_phred to null.
    # we do this so that a single group's output file won't be huge/won't fit in memory if one or a handful of paths
    # have abnormally high number of corresponding reads
    idx = pl.int_range(pl.len(), dtype=pl.UInt64)
    return df.with_columns(
        **{
            c: pl.when(idx < size)
            .then(pl.col(c))
            .otherwise(pl.lit(fill_value))
            .over(group_by)
            # these are the columns containing most of the bytes
            for c in columns
        }
    )
    # if we wanted to entirely filter out rows beyond limit_depth in each grouping, we could do this instead:
    # df = df.select(pl.all().head(limit_depth).over(hash_expr, mapping_strategy="explode"))


def compute_consensus_seqs(
    input_filename,
    output_filename=None,
    fasta_filename=None,
    fastq_filename=None,
    input_format=None,
    output_format=None,
    group=None,
    hash_column=None,
    min_depth=None,
    min_duplex_depth=None,
    limit_depth=None,
    max_length=None,
    method="abpoa",
    use_phreds=None,
    output_phreds=None,
    consensus_kwargs={},
    include=None,
    include_prefix=None,
    exclude=None,
    exclude_prefix=None,
    max_divergence=None,
    skip_consensus=False,
    flag_columns=DEFAULT_FLAG_COLUMNS,
    segments_struct="grouping_segments",
    variant_sep="=",
):
    if not input_filename:
        return  # no op if given no input
    if fastq_filename and not output_phreds:
        raise ValueError("cannot output fastq without --output-phreds")
    # detect format based on first of potentially many filenames/glob patterns
    input_format = detect_format(
        input_format, input_filename[0], ["arrow", "parquet"], glob=True
    )
    if output_filename:
        output_format = detect_format(
            output_format, output_filename, ["arrow", "parquet"]
        )
    if not any([output_filename, fasta_filename, fastq_filename]):
        raise ValueError("at least one of --output, --fasta, --fastq must be given")
    with pl.StringCache():
        # TODO: waiting on Polars to support streaming this query
        if input_format == "arrow":
            df_inputs = [pl.scan_ipc(f) for f in input_filename]
        elif input_format == "parquet":
            df_inputs = [pl.scan_parquet(f) for f in input_filename]
        df_columns = df_inputs[0].collect_schema().names()
        if hash_column is not None and hash_column in df_columns:
            hash_expr = pl.col(hash_column)
        else:
            hash_expr = categorical_list_hash(pl.col("path"))
        ###
        # TODO: memory management issue in polars/jemalloc: https://github.com/pola-rs/polars/issues/19497
        # so instead of executing the entire query lazily, we prefilter all input tables
        # one-by-one for rows with the correct group, then concat those dataframes
        # and execute the rest of the query on that
        # once this is fixed, the following line can replace everything between BEGIN/END WORKAROUND
        # df = pl.concat(df_inputs, how="diagonal")
        ### BEGIN WORKAROUND
        dfs = []
        for df in df_inputs:
            if group:
                df = df.filter(hash_expr % group[1] == group[0])
            if limit_depth:
                # if there's a single group with a large number of rows
                # we do this for each input file individually so that we bound memory usage above by
                # limit_depth * num_input_files * mean_size_per_row
                # we call limit_group_content again below to ensure that the final output
                # is bounded by limit_depth * mean_size_per_row
                # TODO: note that if a single (too-large-for-memory) input file has a large number of rows
                # that all belong to a single oversized group, we may still run into memory issues
                # this should be fixed by using df.collect(streaming=True) below once
                # the streaming engine is implemented for this query
                # we filter first so that we exclude all rows that were already nulled-out
                # during loading each input file above
                df = df.filter(*[pl.col(c).is_not_null() for c in LIMIT_DEPTH_COLUMNS])
                df = limit_group_content(
                    df, LIMIT_DEPTH_COLUMNS, hash_expr, limit_depth
                )
            if max_length:
                df = df.filter(pl.col("read_seq").str.len_bytes() <= max_length)
            # TODO: round-tripping through arrow
            df = pl.from_arrow(df.collect().to_arrow()).lazy()
            # when bug is fixed, we should do this instead:
            # df = df.collect().lazy()
            dfs.append(df)
        df = pl.concat(dfs, how="diagonal")
        ### END WORKAROUND
        if not skip_consensus:
            for col in flag_columns:
                if col in df_columns:
                    df = df.filter(pl.col(col))
        # redundant now, but uncomment when above length filtering is removed
        # if max_length:
        #     df = df.filter(pl.col("read_seq").str.len_bytes() <= max_length)
        if max_divergence is not None:
            if "max_divergence" not in df_columns:
                all_segments = unique_segments(df, pl.col("path"))
                selected_segments = filter_segments(
                    all_segments,
                    include=include,
                    include_prefix=include_prefix,
                    exclude=exclude,
                    exclude_prefix=exclude_prefix,
                )
                if variant_sep:
                    selected_segments = list(
                        set([s.split(variant_sep)[0] for s in selected_segments])
                    )
                df = compute_divergences(
                    df, selected_segments, struct_name=segments_struct
                )
            df = df.filter(pl.col("max_divergence") <= max_divergence)
        # we want to filter (--min-depth, --min-duplex-depth) on the number of reads
        # that may be used to construct consensuses (up to --limit-depth), so we need
        # to compute grouping_depth after filtering on --max-length, --max-divergence
        # in principle we could also record depths before filtering, but this is not terribly
        # useful and you can just go back to the output of PREPARE_CONSENSUS and look up group
        # size for a given path_hash
        depth_columns = compute_depth(df, over="path", prefix="grouping_")
        if not skip_consensus:
            df = df.with_columns(**depth_columns)
        if min_depth:
            df = df.filter(depth_columns["grouping_depth"] >= min_depth)
        if min_duplex_depth and "grouping_duplex_depth" in depth_columns:
            df = df.filter(depth_columns["grouping_duplex_depth"] >= min_duplex_depth)
        if skip_consensus and limit_depth:
            df = limit_group_content(df, LIMIT_DEPTH_COLUMNS, hash_expr, limit_depth)
        column_order = [
            "name",
            "path",
            "path_hash",
            "read_seq",
            "read_phred",
            "reverse_complement",
            "rq",
            "qs",
            "dx",
            "grouping_depth",
            "grouping_duplex_depth",
            segments_struct,
        ]
        if skip_consensus and flag_columns:
            column_order.extend(flag_columns)
        columns = set(column_order)
        if not skip_consensus:
            columns.discard(segments_struct)
        df_columns = df.collect_schema().names()
        # set intersections don't preserve order
        columns &= set(df_columns)
        if "rq" in columns:
            # PacBio CCS uses the "qs" tag for something else, so ignore if "rq" is present
            columns.discard("qs")
        df = df.select(pl.col(sorted(columns, key=column_order.index)))
        if not skip_consensus:
            df = df.filter(pl.col("read_seq").is_not_null())
            # path is already included by group_by
            agg_columns = sorted(
                set(["path_hash", "grouping_depth", "grouping_duplex_depth"]) & columns,
                key=column_order.index,
            )
            df = map_read_groups(
                df,
                partial(
                    get_consensus_group_by,
                    method=method,
                    use_phreds=use_phreds,
                    return_phreds=output_phreds,
                    **consensus_kwargs,
                ),
                pl.col(agg_columns).first(),
                *compute_depth(df, prefix="consensus_").values(),
                max_group_size=limit_depth,
            )
        # TODO: try streaming?
        df = df.collect()
        if not skip_consensus:
            if fasta_filename:
                write_fastx(
                    fasta_filename,
                    df.get_column("consensus_seq"),
                    names=df.get_column("name"),
                )
            if fastq_filename:
                write_fastx(
                    fastq_filename,
                    df.get_column("consensus_seq"),
                    phreds=df.get_column("consensus_phred"),
                    names=df.get_column("name"),
                )
        if output_format == "arrow":
            df.write_ipc(output_filename)
            # TODO: try streaming?
            # df.sink_ipc(output_filename)
        elif output_format == "parquet":
            df.write_parquet(output_filename)
            # TODO: try streaming?
            # df.sink_parquet(output_filename)


def _parse_group(ctx, param, value):
    if value:
        try:
            assert value.count("/") == 1
            return tuple(int(num.strip()) for num in value.split("/"))
        except (AssertionError, ValueError):
            raise click.BadParameter("expecting --group <group_id>/<num_groups>")
    else:
        return None


@click.command(context_settings={"show_default": True})
@click.option("--output", type=click.Path())
@click.option("--fasta", type=click.Path())
@click.option("--fastq", type=click.Path())
@click.option(
    "-i",
    "--input-format",
    type=click.Choice(["parquet", "arrow"], case_sensitive=False),
)
@click.option(
    "-o",
    "--output-format",
    type=click.Choice(["parquet", "arrow"], case_sensitive=False),
)
@click.option("--group", type=str, callback=_parse_group)
@click.option("--hash-col", default="path_hash")
@click.option("--no-hash-col", is_flag=True)
@click.option("--min-depth", type=int)
@click.option("--min-duplex-depth", type=int)
@click.option("--limit-depth", type=int, default=50)  # TODO: is this a good default?
@click.option("--no-limit_depth", is_flag=True)
@click.option("--max-length", type=int, default=20_000)
@click.option("--no-max-length", is_flag=True)
@click.option(
    "--method", type=click.Choice(["abpoa", "spoa", "first"]), default="abpoa"
)
@click.option("--phred-input/--no-phred-input", default=False)  # TODO
@click.option("--phred-output/--no-phred-output", default=True)  # TODO
@click.option("-p", "--param", type=(str, str), multiple=True, callback=parse_kv)
@filter_gfa_options
@click.option("--max-divergence", type=float)
@click.option("--skip-consensus", is_flag=True)
@click.option(
    "--filter-flag",
    default=DEFAULT_FLAG_COLUMNS,
    multiple=True,
    callback=split_delimited_list,
)
@click.option("--no-filter-flag", is_flag=True)
@click.option("--segments-struct", default="grouping_segments")
@click.option("--variant-sep", default="=")
@click.option("--no-variant-sep", is_flag=True)
@click.argument("input", type=str, nargs=-1)
def cli(
    input,
    output,
    fasta,
    fastq,
    input_format,
    output_format,
    group,
    hash_col,
    no_hash_col,
    min_depth,
    min_duplex_depth,
    limit_depth,
    no_limit_depth,
    max_length,
    no_max_length,
    method,
    phred_input,
    phred_output,
    param,
    include,
    include_prefix,
    exclude,
    exclude_prefix,
    max_divergence,
    skip_consensus,
    filter_flag,
    no_filter_flag,
    segments_struct,
    variant_sep,
    no_variant_sep,
):
    compute_consensus_seqs(
        input,
        output,
        fasta,
        fastq,
        input_format,
        output_format,
        group,
        None if no_hash_col else hash_col,
        min_depth,
        min_duplex_depth,
        None if no_limit_depth else limit_depth,
        None if no_max_length else max_length,
        method,
        phred_input,
        phred_output,
        param,
        include,
        include_prefix,
        exclude,
        exclude_prefix,
        max_divergence,
        skip_consensus,
        [] if no_filter_flag else filter_flag,
        segments_struct,
        None if no_variant_sep else variant_sep,
    )


if __name__ == "__main__":
    cli()
