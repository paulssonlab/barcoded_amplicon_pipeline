#!/usr/bin/env python
import sys
from functools import partial
from pathlib import Path

import click
import polars as pl

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.consensus import get_consensus_group_by
from paulssonlab.sequencing.io import write_fastx
from paulssonlab.sequencing.processing import compute_depth, map_read_groups
from paulssonlab.sequencing.util import detect_format


def compute_consensus_seqs(
    input_filename,
    output_filename=None,
    fasta_filename=None,
    fastq_filename=None,
    input_format=None,
    output_format=None,
    group=None,
    min_depth=None,
    min_simplex_depth=None,
    min_duplex_depth=None,
):
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
    if input_format == "arrow":
        df = pl.concat([pl.scan_ipc(f) for f in input_filename])
    elif input_format == "parquet":
        df = pl.concat([pl.scan_parquet(f) for f in input_filename])
    if group:
        df = df.filter(pl.col("path").hash() % group[1] == group[0])
    df = compute_depth(df)
    exprs = []
    if min_depth:
        exprs.append(pl.col("depth") > min_depth)
    if min_simplex_depth:
        exprs.append(pl.col("simplex_depth") > min_simplex_depth)
    if min_duplex_depth:
        exprs.append(pl.col("duplex_depth") > min_duplex_depth)
    if exprs:
        df = df.filter(*exprs)
    df = map_read_groups(df, partial(get_consensus_group_by, return_phreds=True))
    # TODO: try streaming?
    df = df.collect()
    if fasta_filename:
        write_fastx(fasta_filename, df.get_column("consensus_seq"))
    if fastq_filename:
        write_fastx(fastq_filename, df.get_column("consensus_phred"))
    if output_format == "arrow":
        df.write_ipc(output_filename)
    elif output_format == "parquet":
        df.write_parquet(output_filename)


def _parse_group(ctx, param, value):
    try:
        assert value.count("/") == 1
        return tuple(int(num.strip()) for num in value.split("/"))
    except:
        raise click.BadParameter("expecting --group <group_id>/<num_groups>")


@click.command()
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
@click.option("--min-depth", type=int)
@click.option("--min-simplex-depth", type=int)
@click.option("--min-duplex-depth", type=int)
@click.argument("input", type=str, nargs=-1)
def cli(
    input,
    output,
    fasta,
    fastq,
    input_format,
    output_format,
    group,
    min_depth,
    min_simplex_depth,
    min_duplex_depth,
):
    compute_consensus_seqs(
        input,
        output,
        fasta,
        fastq,
        input_format,
        output_format,
        group,
        min_depth,
        min_simplex_depth,
        min_duplex_depth,
    )


if __name__ == "__main__":
    cli()
