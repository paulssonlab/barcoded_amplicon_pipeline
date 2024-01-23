#!/usr/bin/env python
import sys
from pathlib import Path

import click
import gfapy

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.gfa import filter_gfa, filter_gfa_options


@click.command()
@filter_gfa_options
@click.option("--mapping", type=click.Path())
@click.argument("input", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path())
def cli(input, output, mapping, include, include_prefix, exclude, exclude_prefix):
    gfa = gfapy.Gfa.from_file(input)
    gfa = filter_gfa(gfa, include, include_prefix, exclude, exclude_prefix)
    name_mapping = []
    for idx, segment in enumerate(gfa.segments):
        name_mapping.append((idx, segment.name))
        segment.name = str(idx)
    gfa.to_file(output)
    if mapping:
        with open(mapping, "w") as f:
            for idx, name in name_mapping:
                f.write(f"{idx}\t{name}\n")
            f.write("\n")


if __name__ == "__main__":
    cli()
