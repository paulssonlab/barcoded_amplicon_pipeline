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
from paulssonlab.sequencing.processing import identify_usable_reads, normalize_paths
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
    keep_partial_paths,
    hash_paths,
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
    if keep_partial_paths:
        endpoints = None
    else:
        endpoints = dag_endpoints(graph, wccs=wccs)
    with pl.StringCache():
        if input_format == "arrow":
            df = pl.concat([pl.scan_ipc(f) for f in input_filename])
        elif input_format == "parquet":
            df = pl.concat([pl.scan_parquet(f) for f in input_filename])
        df = normalize_paths(
            df, forward_segments, endpoints=endpoints, hash_paths=hash_paths
        )
        # identify_usable_reads is much faster when working on in-memory data
        df = df.collect().lazy()
        # there's no difference in doing this lazily or not,
        # so I'm arbitrarily choosing to do it lazily
        df = identify_usable_reads(df)
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
@click.option(
    "--keep-partial-paths/--no-keep-partial-paths",
    default=False,
    help="Keep alignments with paths that do not span graph end-to-end",
)
@click.option(
    "--hash-paths/--no-hash-paths", default=True, help="Precompute path hashes"
)
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
    keep_partial_paths,
    hash_paths,
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
        keep_partial_paths,
        hash_paths,
    )


if __name__ == "__main__":
    cli()
