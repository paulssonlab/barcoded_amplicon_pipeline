#!/usr/bin/env python
import sys
from pathlib import Path
from urllib.parse import quote_plus

import click
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.csv as csv

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.util.cli import split_delimited_list


def write_read_ids(tsv_filename, output_dir, fields, chunks):
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
                    if chunks != 0:
                        key = hash(key) % chunks
                        key_str = f"chunk={key}"
                    if not (output := output_files.get(key)):
                        if chunks == 0:
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


@click.command()
@click.option(
    "-F",
    "--fields",
    default=[],
    multiple=True,
    show_default=True,
    callback=split_delimited_list,
)
@click.option(
    "-c", "--chunks", type=click.IntRange(min=0), default=0, show_default=True
)
@click.argument("input_tsv", type=click.Path(exists=True, dir_okay=False))
@click.argument("output_dir", type=click.Path())
def cli(input_tsv, output_dir, fields, chunks):
    write_read_ids(input_tsv, output_dir, fields, chunks)


if __name__ == "__main__":
    cli()
