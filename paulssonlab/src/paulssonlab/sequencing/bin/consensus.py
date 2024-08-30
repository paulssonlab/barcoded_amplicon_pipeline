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
            df_input = pl.concat(
                [pl.scan_ipc(f) for f in input_filename], how="diagonal"
            )
        elif input_format == "parquet":
            df_input = pl.concat(
                [pl.scan_parquet(f) for f in input_filename], how="diagonal"
            )
        df_input = df_input.with_columns(_idx=pl.int_range(pl.len(), dtype=pl.UInt64))
        df = df_input
        # TODO
        # df_columns = df.collect_schema().names()
        df_columns = df.columns
        if not skip_consensus:
            for col in flag_columns:
                if col in df_columns:
                    df = df.filter(pl.col(col))
        if hash_column is not None and hash_column in df_columns:
            hash_expr = pl.col(hash_column)
        else:
            hash_expr = categorical_list_hash(pl.col("path"))
        if group:
            df = df.filter(hash_expr % group[1] == group[0])
        if max_length:
            df = df.filter(pl.col("read_seq").str.len_bytes() <= max_length)
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
        depth_columns = compute_depth(df, over="path", prefix="grouping_")
        if not skip_consensus:
            df = df.with_columns(**depth_columns)
        if min_depth:
            df = df.filter(depth_columns["grouping_depth"] >= min_depth)
        if min_duplex_depth and "grouping_duplex_depth" in depth_columns:
            df = df.filter(depth_columns["grouping_duplex_depth"] >= min_duplex_depth)
        if skip_consensus and limit_depth:
            # if we are not computing consensuses, we still include more than limit_depth number of reads
            # but set their read_seq and read_phred to null.
            # we do this so that a single group's output file won't be huge/won't fit in memory if one or a handful of paths
            # have abnormally high number of corresponding reads
            # TODO
            # the obvious way to do this is:
            # df = df.select(pl.all().head(limit_depth).over(hash_expr, mapping_strategy="explode"))
            # but this doesn't reduce memory usage. I also tried:
            # idx = pl.int_range(pl.len(), dtype=pl.UInt32)
            # df = df.with_columns(
            #     **{c: pl.when(idx < limit_depth).then(pl.col(c)).otherwise(pl.lit(None)).over(hash_expr, mapping_strategy="explode") for c in ("grouping_segments", "read_seq", "read_phred")}
            # )
            # which also did not help. so instead we do this (see also .join(...) below)
            df = df.select(
                pl.all()
                .exclude(segments_struct, "read_seq", "read_phred")
                .head(limit_depth)
                .over(hash_expr, mapping_strategy="explode")
            )
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
        # TODO
        # df_columns = df.collect_schema().names()
        df_columns = df.columns
        # set intersections don't preserve order
        columns &= set(df_columns)
        if "rq" in columns:
            # PacBio CCS uses the "qs" tag for something else, so ignore if "rq" is present
            columns.discard("qs")
        if skip_consensus and limit_depth:
            df_join = df_input.select(
                pl.col("_idx", segments_struct, "read_seq", "read_phred")
            )  # .head(1_000).collect()
            df = df.join(df_join, on="_idx", how="left", coalesce=True)  # .lazy()
            df = df.select(pl.col(sorted(columns, key=column_order.index)))
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
@click.option("--limit-depth", type=int, default=50)  # TODO?
@click.option("--max-length", type=int)
@click.option("--method", type=click.Choice(["abpoa", "spoa"]), default="abpoa")
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
    max_length,
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
        limit_depth,
        max_length,
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
