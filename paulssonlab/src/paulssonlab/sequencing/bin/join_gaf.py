#!/usr/bin/env python
import sys
from pathlib import Path

import click
import polars as pl
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.io import iter_bam_and_gaf, read_gaf
from paulssonlab.sequencing.processing import join_dfs
from paulssonlab.sequencing.util import detect_format


# TODO: obsolete when dorado 0.4.1 is released
def fix_dx(table):
    if "dx" in table.column_names:
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
    input_filename,
    output_filename,
    input_format,
    output_format,
    gaf_filename,
    include_unaligned,
    reads_prefix,
    gaf_prefix,
):
    input_format = detect_format(
        input_format,
        input_filename[0],
        ["bam", "arrow", "parquet"],
        glob=True,
    )
    output_format = detect_format(output_format, output_filename, ["arrow", "parquet"])
    if input_format == "bam":
        if len(input_filename) != 1:
            raise ValueError("joining more than one bam file is unsupported")
        batches = iter_bam_and_gaf(
            input_filename[0], gaf_filename, include_unaligned=include_unaligned
        )
        table = pa.Table.from_batches(batches)
        # can't do this streaming because we need to unify dictionaries
        table = table.unify_dictionaries()
    else:
        if input_format == "arrow":
            reads = pl.concat([pl.scan_ipc(f) for f in input_filename])
        elif input_format == "parquet":
            reads = pl.concat([pl.scan_parquet(f) for f in input_filename])
        gaf = pl.from_arrow(read_gaf(gaf_filename)).lazy()
        df = join_dfs(
            reads,
            gaf,
            on="name",
            how="left",
            left_prefix=reads_prefix,
            right_prefix=gaf_prefix,
        )
        table = df.collect().to_arrow()
    table = fix_dx(table)  # TODO
    if output_format == "arrow":
        with pa.ipc.new_file(
            output_filename,
            table.schema,
        ) as writer:
            writer.write(table)
    elif output_format == "parquet":
        pq.write_table(table, output_filename)
    else:
        raise ValueError(f"invalid format: {output_format}")


@click.command(context_settings={"show_default": True})
@click.option("--gaf", type=click.Path(exists=True, dir_okay=False), required=True)
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
@click.option("--reads-prefix", type=str)
@click.option("--gaf-prefix", type=str)
@click.argument("input", type=str, nargs=-1)
@click.argument("output", type=click.Path())
def cli(
    input,
    output,
    input_format,
    output_format,
    gaf,
    include_unaligned,
    reads_prefix,
    gaf_prefix,
):
    join_reads_and_gaf(
        input,
        output,
        input_format,
        output_format,
        gaf,
        include_unaligned,
        reads_prefix,
        gaf_prefix,
    )


if __name__ == "__main__":
    cli()
