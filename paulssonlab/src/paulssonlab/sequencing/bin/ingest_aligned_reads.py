#!/usr/bin/env python
import sys
from pathlib import Path

import click
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
    table = fix_dx(table)  # TODO
    if format == "arrow":
        with pa.ipc.new_file(
            output_filename,
            table.schema,
            options=pa.ipc.IpcWriteOptions(compression="lz4"),
        ) as writer:
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
