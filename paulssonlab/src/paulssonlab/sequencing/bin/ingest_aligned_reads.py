#!/usr/bin/env python
import sys
from pathlib import Path

import click
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.io import iter_bam_and_gaf


def write_bam_and_gaf(
    bam_filename, gaf_filename, output_filename, format, include_unaligned
):
    if format is None:
        if output_filename.endswith(".arrow"):
            format = "arrow"
        else:
            format = "parquet"
    table = pa.Table.from_batches(iter_bam_and_gaf(bam_filename, gaf_filename))
    # can't do this streaming because we need to unify dictionaries
    table = table.unify_dictionaries()
    if format == "arrow":
        with pa.ipc.new_file(output_filename, table.schema) as writer:
            writer.write(table)
    elif format == "parquet":
        pq.write_table(table, output_filename)
    else:
        raise ValueError(f"unknown format: {format}")


@click.command()
@click.option(
    "-f",
    "--format",
    type=click.Choice(["parquet", "arrow"], case_sensitive=False),
)
@click.option("--include-unaligned/--no-include-unaligned", default=True)
@click.argument("bam", type=click.Path(exists=True, dir_okay=False))
@click.argument("gaf", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path(exists=False))
def cli(bam, gaf, output, format, include_unaligned):
    write_bam_and_gaf(bam, gaf, output, format, include_unaligned)


if __name__ == "__main__":
    cli()
