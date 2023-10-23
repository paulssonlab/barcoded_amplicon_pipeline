#!/usr/bin/env python
import sys
from pathlib import Path

import click
import gfapy

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.gfa import filter_gfa, filter_gfa_options


@click.command()
@filter_gfa_options
@click.argument("input", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path())
def cli(input, output, include, include_prefix, exclude, exclude_prefix):
    gfa = gfapy.Gfa.from_file(input)
    gfa = filter_gfa(gfa, include, include_prefix, exclude, exclude_prefix)
    gfa.to_file(output)


if __name__ == "__main__":
    cli()
