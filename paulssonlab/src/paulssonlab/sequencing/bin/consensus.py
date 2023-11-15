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
from paulssonlab.util.cli import parse_kv


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
    min_simplex_depth=None,
    min_duplex_depth=None,
    limit_depth=None,
    method="abpoa",
    use_phreds=None,
    output_phreds=None,
    consensus_kwargs={},
    skip_consensus=False,
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
    # TODO
    import time

    start = time.time()
    with pl.StringCache():
        # TODO: waiting on Polars to support streaming this query
        if input_format == "arrow":
            df = pl.concat([pl.scan_ipc(f) for f in input_filename])
        elif input_format == "parquet":
            df = pl.concat([pl.scan_parquet(f) for f in input_filename])
        if group:
            if hash_column is not None and hash_column in df.columns:
                hash_expr = pl.col(hash_column)
            else:
                hash_expr = pl.col("path").hash()
            df = df.filter(hash_expr % group[1] == group[0])
        df = df.select(
            pl.col(
                "name",
                "path",
                "is_duplex",
                "read_seq",
                "read_phred",
                "reverse_complement",
            )
        )
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
        if not skip_consensus:
            df = map_read_groups(
                df,
                partial(
                    get_consensus_group_by,
                    method=method,
                    use_phreds=use_phreds,
                    return_phreds=output_phreds,
                    **consensus_kwargs,
                ),
                max_group_size=limit_depth,
            )
        # TODO: try streaming?
        df = df.collect()
        print(f"GROUP BY DONE {time.time() - start}")
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
        print(f"FASTX WRITE DONE {time.time() - start}")
        if output_format == "arrow":
            df.write_ipc(output_filename)
        elif output_format == "parquet":
            df.write_parquet(output_filename)
        print(f"TABULAR WRITE DONE {time.time() - start}")


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
@click.option("--min-simplex-depth", type=int)
@click.option("--min-duplex-depth", type=int)
@click.option("--limit-depth", type=int, default=50)
@click.option("--method", type=click.Choice(["abpoa", "spoa"]), default="abpoa")
@click.option("--phred-input/--no-phred-input", default=False)  # TODO
@click.option("--phred-output/--no-phred-output", default=True)  # TODO
@click.option("-p", "--param", type=(str, str), multiple=True, callback=parse_kv)
@click.option("--skip-consensus", is_flag=True)
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
    min_simplex_depth,
    min_duplex_depth,
    limit_depth,
    method,
    phred_input,
    phred_output,
    param,
    skip_consensus,
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
        min_simplex_depth,
        min_duplex_depth,
        limit_depth,
        method,
        phred_input,
        phred_output,
        param,
        skip_consensus,
    )


if __name__ == "__main__":
    cli()
