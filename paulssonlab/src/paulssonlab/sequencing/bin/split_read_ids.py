#!/usr/bin/env python
from pathlib import Path
from urllib.parse import quote_plus

import click
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as csv


def write_read_ids(tsv_filename, output_dir, fields):
    if not fields:
        raise ValueError("at least one field must be specified")
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    output_files = {}
    try:
        with csv.open_csv(
            tsv_filename,
            parse_options=csv.ParseOptions(delimiter="\t"),
            read_options=csv.ReadOptions(block_size=100_000_000),
        ) as f:
            while True:
                try:
                    table = pa.Table.from_batches([f.read_next_batch()])
                except StopIteration:
                    break
                # FROM: https://github.com/apache/arrow/issues/14882#issuecomment-1347260487
                if not fields:
                    fields = list(set(table.column_names) - set(["read_id"]))
                    if not fields:
                        raise ValueError("need at least one column to split on")
                groups = table.group_by(fields).aggregate([("read_id", "list")])
                field_values = {field: groups.column(field) for field in fields}
                read_ids = groups.column("read_id_list")
                for idx in range(len(groups)):
                    key = tuple([values[idx] for values in field_values.values()])
                    if not (output := output_files.get(key)):
                        key_str = ",".join(
                            [
                                f"{field}={quote_plus(str(value))}"
                                for field, value in zip(fields, key)
                            ]
                        )
                        output = output_files[key] = open(
                            output_dir / f"{key_str}.txt", "w"
                        )
                    output.write(pc.binary_join(read_ids[idx], "\n").as_py() + "\n")
    finally:
        for file in output_files.values():
            file.close()


def _split_delimited_list(ctx, param, value):
    return [c.strip() for c in value.split(",")]


@click.command()
@click.option(
    "-F",
    "--fields",
    default=[],
    show_default=True,
    callback=_split_delimited_list,
)
@click.argument("input_tsv", type=click.Path(exists=True, dir_okay=False))
@click.argument("output_dir", type=click.Path(exists=False))
def cli(input_tsv, output_dir, fields):
    write_read_ids(input_tsv, output_dir, fields)


if __name__ == "__main__":
    cli()
