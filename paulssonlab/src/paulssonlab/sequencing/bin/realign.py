#!/usr/bin/env python
import sys
from pathlib import Path

import click
import polars as pl
from gfapy import Gfa

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.processing import pairwise_align_to_path
from paulssonlab.sequencing.util import detect_format
from paulssonlab.util.cli import parse_kv


def realign(
    gfa_filename,
    input_filename,
    output_filename,
    input_format,
    output_format,
    align_kwargs={},
):
    if not input_filename:
        return  # no op if given no input
    # detect format based on first of potentially many filenames/glob patterns
    input_format = detect_format(
        input_format, input_filename[0], ["arrow", "parquet"], glob=True
    )
    if output_filename:
        output_format = detect_format(
            output_format, output_filename, ["arrow", "parquet"]
        )
    gfa = Gfa.from_file(gfa_filename)
    if input_format == "arrow":
        df = pl.concat([pl.scan_ipc(f) for f in input_filename])
    elif input_format == "parquet":
        df = pl.concat([pl.scan_parquet(f) for f in input_filename])
    df = df.collect()
    df = pairwise_align_to_path(df, gfa, **align_kwargs)
    if output_format == "arrow":
        df.write_ipc(output_filename)
    elif output_format == "parquet":
        df.write_parquet(output_filename)


def _parse_params(ctx, param, value):
    d = {}
    for k, v in value:
        try:
            v = int(v)
        except:
            try:
                v = float(v)
            except:
                pass
        d[k] = v
    return d


@click.command()
@click.option("--gfa", type=click.Path(exists=True, dir_okay=False), required=True)
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
@click.option("--method", type=click.Choice(["parasail", "pywfa"]), default="parasail")
@click.option("--degenerate/--no-degenerate", default=True)
@click.option("-p", "--param", type=(str, str), multiple=True, callback=parse_kv)
@click.argument("input", type=click.Path(exists=True, dir_okay=False), nargs=-1)
@click.argument("output", type=click.Path())
def cli(gfa, input, output, input_format, output_format, method, degenerate, param):
    realign(
        gfa,
        input,
        output,
        input_format,
        output_format,
        {"method": method, "degenerate": degenerate, **param},
    )


if __name__ == "__main__":
    cli()
