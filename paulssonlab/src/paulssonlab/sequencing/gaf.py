import re

import pyarrow as pa
import pyarrow.compute as pc
from pyarrow import csv

# SEE: http://samtools.github.io/hts-specs/SAMv1.pdf
# and https://samtools.github.io/hts-specs/SAMtags.pdf
# pyarrow CSV parser only supports pa.dictionary with int32 indices
SAM_TAG_TYPES = {
    "A": pa.dictionary(pa.int32(), pa.string()),
    "f": pa.float32(),
    "i": pa.int32(),
    "Z": pa.string(),
}
GAF_COLUMN_TYPES = {
    "query_length": pa.uint64(),
    "query_start": pa.uint64(),
    "query_end": pa.uint64(),
    "strand": pa.dictionary(pa.int32(), pa.string()),
    "path": pa.string(),
    "path_length": pa.uint64(),
    "path_start": pa.uint64(),
    "path_end": pa.uint64(),
    "residue_matches": pa.uint64(),
    "block_length": pa.uint64(),
    "mapping_quality": pa.uint8(),
}
SAM_TAG_REGEX = re.compile(
    r"^(?P<tag>[a-zA-Z0-9]+):(?P<tag_value>A:.|f:\d+(\.\d+)?|i:\d+|Z:.*)$"
)


def parse_gaf_types(gaf_filename):
    with open(gaf_filename, "r") as f:
        first_row = f.readline().split("\t")
    columns_to_parse = {}
    column_types = []
    for idx in reversed(range(len(first_row))):
        if match := SAM_TAG_REGEX.match(first_row[idx]):
            tag = match.group("tag")
            column_types.append((tag, pa.string()))
            tag_value = match.group("tag_value")
            columns_to_parse[tag] = tag_value[: tag_value.index(":")]
        else:
            break
    column_types.extend(reversed(GAF_COLUMN_TYPES.items()))
    for idx in reversed(range(idx + 1 - len(GAF_COLUMN_TYPES))):
        if match := SAM_TAG_REGEX.match(first_row[idx]):
            tag = match.group("tag")
            column_types.append((tag, pa.string()))
            tag_value = match.group("tag_value")
            type_ = tag_value[: tag_value.index(":")]
            columns_to_parse[tag] = type_
        else:
            if idx != 0:
                raise ValueError("expecting SAM tags following FASTQ read name")
            else:
                column_types.append(("name", pa.string()))
    column_types = dict(reversed(column_types))
    return column_types, columns_to_parse


def parse_gaf_table(table, columns_to_parse):
    # TODO: we could convert string read UUIDs (and semicolon-delimited pairs of UUIDs)
    # to an extension type to save a small amount of space
    # SEE: https://arrow.apache.org/docs/python/extending_types.html#defining-extension-types-user-defined-types
    for tag, type_ in columns_to_parse.items():
        col_idx = table.column_names.index(tag)
        new_column = pc.replace_substring_regex(table[tag], f"{tag}:{type_}:", "").cast(
            SAM_TAG_TYPES[type_]
        )
        table = table.set_column(col_idx, tag, new_column)
    path = pa.array(
        [re.split(r"(?=<|>)", s.as_py())[1:] for s in table.column("path")],
        type=pa.list_(pa.dictionary(pa.int16(), pa.string())),
    )
    table = table.set_column(table.column_names.index("path"), "path", path)
    return table


def iter_gaf(gaf_filename, block_size=None):
    column_types, columns_to_parse = parse_gaf_types(gaf_filename)
    read_options = csv.ReadOptions(
        column_names=column_types.keys(), block_size=block_size
    )
    parse_options = csv.ParseOptions(delimiter="\t")
    convert_options = csv.ConvertOptions(column_types=column_types)
    with csv.open_csv(
        gaf_filename,
        read_options=read_options,
        parse_options=parse_options,
        convert_options=convert_options,
    ) as f:
        while True:
            try:
                table = parse_gaf_table(
                    pa.Table.from_batches([f.read_next_batch()]), columns_to_parse
                )
            except StopIteration:
                break
            yield table


def read_gaf(gaf_filename, **kwargs):
    return pa.concat_tables(list(iter_gaf(gaf_filename, **kwargs)))
