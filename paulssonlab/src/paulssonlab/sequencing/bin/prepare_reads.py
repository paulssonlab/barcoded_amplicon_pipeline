#!/usr/bin/env python
import sys
from pathlib import Path

import click
import gfapy
import networkx as nx
import polars as pl

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.gfa import (
    dag_endpoints,
    dag_forward_segments,
    filter_gfa,
    filter_gfa_options,
    gfa_to_dag,
)
from paulssonlab.sequencing.processing import prepare_reads as _prepare_reads
from paulssonlab.sequencing.util import detect_format


def prepare_reads(
    gfa_filename,
    input_filename,
    output_filename,
    input_format,
    output_format,
    include,
    include_prefix,
    exclude,
    exclude_prefix,
):
    input_format = detect_format(
        input_format,
        input_filename[0],
        ["arrow", "parquet"],
        glob=True,
    )
    output_format = detect_format(output_format, output_filename, ["arrow", "parquet"])
    gfa = gfapy.Gfa.from_file(gfa_filename)
    gfa = filter_gfa(gfa, include, include_prefix, exclude, exclude_prefix)
    graph = gfa_to_dag(gfa)
    # weakly_connected_components is a generator, so only compute once
    wccs = list(nx.weakly_connected_components(graph))
    forward_segments = dag_forward_segments(graph, wccs=wccs)
    endpoints = dag_endpoints(graph, wccs=wccs)
    with pl.StringCache():
        if input_format == "arrow":
            df = pl.concat([pl.scan_ipc(f) for f in input_filename])
        elif input_format == "parquet":
            df = pl.concat([pl.scan_parquet(f) for f in input_filename])
        df = _prepare_reads(df, forward_segments, endpoints)
        df = df.collect()
        if output_format == "arrow":
            df.write_ipc(output_filename)
        elif output_format == "parquet":
            df.write_parquet(output_filename)


@click.command(context_settings={"show_default": True})
@click.option(
    "-i",
    "--input-format",
    type=click.Choice(["bam", "parquet", "arrow"], case_sensitive=False),
)
@click.option(
    "-o",
    "--output-format",
    type=click.Choice(["parquet", "arrow"], case_sensitive=False),
)
@filter_gfa_options
@click.option("--gfa", type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument("input", type=click.Path(exists=True, dir_okay=False), nargs=-1)
@click.argument("output", type=click.Path())
def cli(
    gfa,
    input,
    output,
    input_format,
    output_format,
    include,
    include_prefix,
    exclude,
    exclude_prefix,
):
    prepare_reads(
        gfa,
        input,
        output,
        input_format,
        output_format,
        include,
        include_prefix,
        exclude,
        exclude_prefix,
    )


if __name__ == "__main__":
    cli()
