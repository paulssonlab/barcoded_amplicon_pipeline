#!/usr/bin/env python
import sys
from pathlib import Path

import click
import polars as pl
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.io import iter_bam_and_gaf


# TODO: obsolete when dorado 0.4.1 is released
def fix_dx(table):
    name_col = table.column("name")
    new_dx = pc.if_else(
        pc.is_in(
            name_col,
            pc.list_flatten(
                pc.split_pattern(
                    name_col.filter(pc.match_substring(name_col, ";")), ";"
                )
            ),
        ),
        -1,
        table.column("dx"),
    )
    table = table.set_column(table.column_names.index("dx"), "dx", new_dx)
    return table


def join_reads_and_gaf(
    reads_filename,
    gaf_filename,
    output_filename,
    input_format,
    output_format,
    include_unaligned,
):
    if input_format is None:
        if reads_filename.endswith(".bam"):
            input_format = "bam"
        elif reads_filename.endswith(".arrow"):
            input_format = "arrow"
        elif reads_filename.endswith(".parquet"):
            input_format = "parquet"
        else:
            raise ValueError(f"unknown file extension: {input_filename}")
    if output_format is None:
        if output_filename.endswith(".arrow"):
            format = "arrow"
        elif output_filename.endswith(".parquet"):
            format = "parquet"
        else:
            raise ValueError(f"unknown file extension: {output_filename}")
    if input_format == "bam":
        batches = iter_bam_and_gaf(
            reads_filename, gaf_filename, include_unaligned=include_unaligned
        )
        table = pa.Table.from_batches(batches)
        # can't do this streaming because we need to unify dictionaries
        table = table.unify_dictionaries()
    else:
        if input_format == "arrow":
            reads = pl.scan_ipc(input_filename)
        elif input_format == "parquet":
            reads = pl.scan_parquet(input_filename)
        gaf = pl.from_arrow(read_gaf(gaf_filename))
        df = reads.join(gaf, on="name", how="left")
        table = df.to_arrow()
    table = fix_dx(table)  # TODO
    if format == "arrow":
        with pa.ipc.new_file(
            output_filename,
            table.schema,
        ) as writer:
            writer.write(table)
    elif format == "parquet":
        pq.write_table(table, output_filename)
    else:
        raise ValueError(f"unknown format: {format}")


@click.command()
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
@click.option("--include-unaligned/--no-include-unaligned", default=True)
@click.argument("reads", type=click.Path(exists=True, dir_okay=False))
@click.argument("gaf", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path(exists=False))
def cli(reads, gaf, output, input_format, output_format, include_unaligned):
    join_reads_and_gaf(
        reads,
        gaf,
        output,
        input_format,
        output_format,
        include_unaligned,
    )


if __name__ == "__main__":
    cli()
