#!/usr/bin/env python
import sys
from pathlib import Path

import click
import gfapy
import polars as pl

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.gfa import (
    filter_gfa,
    filter_gfa_options,
    get_forward_segments,
    gfa_to_dag,
)
from paulssonlab.sequencing.processing import identify_usable_reads, normalize_path


def prepare_reads(
    gfa_filename,
    input_filename,
    output_filename,
    format,
    include,
    include_prefix,
    exclude,
    exclude_prefix,
    keep_partial_paths,
):
    if format is None:
        if input_filename.endswith(".arrow"):
            format = "arrow"
        elif input_filename.endswith(".parquet"):
            format = "parquet"
        else:
            raise ValueError(f"unknown file extension: {input_filename}")
    gfa = gfapy.from_file(gfa_filename)
    gfa = filter_gfa()
    with pl.StringCache():
        if format == "arrow":
            df = pl.scan_ipc(input_filename)
        elif format == "parquet":
            df = pl.scan_parquet(input_filename)
        gfa = gfapy.Gfa.from_file(gfa_filename)
        gfa = filter_gfa(gfa, include, include_prefix, exclude, exclude_prefix)
        graph = gfa_to_dag(gfa)
        wccs = nx.weakly_connected_components(graph)
        forward_segments = dag_forward_segments(graph, wccs=wccs)
        if keep_partial_paths:
            endpoints = None
        else:
            endpoints = dag_endpoints(graph, wccs=wccs)
        df = normalize_path(df, forward_segments, endpoints=endpoints)
        # identify_usable_reads is much faster when working on in-memory data
        df = df.collect().lazy()
        df = identify_usable_reads(df)
        df = df.collect()
        if format == "arrow":
            df.write_ipc(output_filename)
        elif format == "parquet":
            df.write_parquet(output_filename)


@click.command()
@click.option(
    "-f",
    "--format",
    type=click.Choice(["parquet", "arrow"], case_sensitive=False),
)
@filter_gfa_options
@click.option(
    "-p",
    "--keep-partial-paths",
    help="Keep alignments with paths that do not span graph end-to-end",
)
@click.argument("gfa", type=click.Path(exists=True, dir_okay=False))
@click.argument("input", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path(exists=False))
def cli(
    gfa,
    input,
    output,
    format,
    include,
    include_prefix,
    exclude,
    exclude_prefix,
    keep_partial_paths,
):
    prepare_reads(
        gfa,
        input,
        output,
        format,
        include,
        include_prefix,
        exclude,
        exclude_prefix,
        keep_partial_paths,
    )


if __name__ == "__main__":
    cli()
