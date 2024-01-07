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
from paulssonlab.sequencing.processing import find_duplex_pairs
from paulssonlab.sequencing.processing import find_duplex_pairs as _find_duplex_pairs
from paulssonlab.sequencing.processing import load_pairing_data


def find_duplex_pairs(
    gfa_filename,
    gaf_filename,
    bam_filename,
    output_filename,
    include,
    include_prefix,
    exclude,
    exclude_prefix,
    end_to_end,
):
    gfa = gfapy.Gfa.from_file(gfa_filename)
    gfa = filter_gfa(gfa, include, include_prefix, exclude, exclude_prefix)
    graph = gfa_to_dag(gfa)
    # weakly_connected_components is a generator, so only compute once
    wccs = list(nx.weakly_connected_components(graph))
    forward_segments = dag_forward_segments(graph, wccs=wccs)
    if end_to_end:
        endpoints = dag_endpoints(graph, wccs=wccs)
    else:
        endpoints = None
    with pl.StringCache():
        df = load_pairing_data(bam_filename, gaf_filename)
        # there's no difference in doing this lazily or not,
        # so I'm arbitrarily choosing to do it lazily
        df = df.lazy()
        df = _find_duplex_pairs(
            df,
            forward_segments,
            endpoints=endpoints,
        )
        df = df.select(["name", "name_right"])
        df = df.collect()
        df.write_csv(output_filename, include_header=False, separator=" ")


@click.command(context_settings={"show_default": True})
@filter_gfa_options
@click.option(
    "--end-to-end/--no-end-to-end",
    default=True,
    help="Only use alignments with paths that span graph end-to-end",
)
@click.option("--gfa", type=click.Path(exists=True, dir_okay=False), required=True)
@click.option("--gaf", type=click.Path(exists=True, dir_okay=False), required=True)
@click.argument("bam", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path())
def cli(
    gfa,
    gaf,
    bam,
    output,
    include,
    include_prefix,
    exclude,
    exclude_prefix,
    end_to_end,
):
    find_duplex_pairs(
        gfa,
        gaf,
        bam,
        output,
        include,
        include_prefix,
        exclude,
        exclude_prefix,
        end_to_end,
    )


if __name__ == "__main__":
    cli()
