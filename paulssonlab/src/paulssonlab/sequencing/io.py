import itertools as it
import re

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pyarrow import csv

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
    r"^(?P<tag>[a-zA-Z0-9]+):(?P<tag_value>A:.|f:(\+|-)?\d+(\.\d+)?|i:(\+|-)?\d+|Z:.*)$"
)
DEFAULT_BAM_TAGS = {
    # ONT
    "RG": pa.dictionary(pa.uint8(), pa.string()),
    "qs": pa.uint8(),
    "ns": pa.int64(),
    "ts": pa.int64(),
    "mx": pa.uint8(),
    "ch": pa.uint16(),
    "rn": pa.uint32(),
    "st": pa.string(),
    "du": pa.float32(),
    "fn": pa.dictionary(pa.uint8(), pa.string()),
    "sm": pa.float32(),
    "sf": pa.float32(),
    "sv": pa.dictionary(pa.uint8(), pa.string()),
    # "mv": ???,
    "dx": pa.int8(),
    "pi": pa.string(),
    "sp": pa.int64(),
    "pt": pa.int64(),
    "MN": pa.int64(),
    # PacBio
}

# SEE: http://samtools.github.io/hts-specs/SAMv1.pdf
# and https://samtools.github.io/hts-specs/SAMtags.pdf
# pyarrow CSV parser only supports pa.dictionary with int32 indices
# B,? (array) types are not used because pysam strips off the subtype
# and only returns "B" as the type
SAM_TAG_TYPES = {
    "A": pa.dictionary(pa.int32(), pa.string()),
    "C": pa.uint8(),  # char (0-255), deprecated?? I'm confused why pysam outputs this
    "B,c": pa.list_(pa.int8()),
    "B,C": pa.list_(pa.uint8()),
    "B,s": pa.list_(pa.int16()),
    "B,S": pa.list_(pa.uint16()),
    "B,i": pa.list_(pa.int32()),
    "B,I": pa.list_(pa.uint32()),
    "B,f": pa.list_(pa.float32()),
    "f": pa.float32(),
    "i": pa.int32(),
    "Z": pa.string(),
}

# python array.typecode to pyarrow types
# I don't understand the difference between h/H and i/I
SAM_TAG_ARRAY_TYPES = {
    "b": pa.list_(pa.int8()),
    "B": pa.list_(pa.uint8()),
    "h": pa.list_(pa.int16()),
    "H": pa.list_(pa.uint16()),
    "i": pa.list_(pa.int16()),
    "I": pa.list_(pa.uint16()),
    "l": pa.list_(pa.int32()),
    "L": pa.list_(pa.uint32()),
    "f": pa.list_(pa.float32()),
}


def pyarrow_type_for_bam(type_, value):
    if type_ in SAM_TAG_TYPES:
        return SAM_TAG_TYPES[type_]
    elif type_ == "B" and value.typecode in SAM_TAG_ARRAY_TYPES:
        return SAM_TAG_ARRAY_TYPES[value.typecode]
    else:
        raise ValueError(f"unknown SAM tag type: {type_} (value: {value})")


def parse_gaf_types(gaf_filename):
    with open(gaf_filename, "r") as f:
        first_row = f.readline().rstrip("\n").split("\t")
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


def parse_gaf_batch(batch, columns_to_parse):
    # TODO: we could convert string read UUIDs (and semicolon-delimited pairs of UUIDs)
    # to an extension type to save a small amount of space
    # SEE: https://arrow.apache.org/docs/python/extending_types.html#defining-extension-types-user-defined-types
    columns = {}
    for col_name, column in zip(batch.column_names, batch.columns):
        if type_ := columns_to_parse.get(col_name):
            columns[col_name] = pc.replace_substring_regex(
                column, f"{col_name}:{type_}:", ""
            ).cast(SAM_TAG_TYPES[type_])
        else:
            columns[col_name] = column
    columns["path"] = pa.array(
        [re.split(r"(?=<|>)", s.as_py())[1:] for s in columns["path"]],
        type=pa.list_(pa.dictionary(pa.int16(), pa.string())),
    )
    return pa.RecordBatch.from_arrays(
        list(columns.values()), names=list(columns.keys())
    )


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
                batch = parse_gaf_batch(f.read_next_batch(), columns_to_parse)
            except StopIteration:
                break
            yield batch


def read_gaf(gaf_filename, **kwargs):
    return pa.Table.from_batches(iter_gaf(gaf_filename, **kwargs))


def iter_bam_and_gaf(
    bam_filename,
    gaf_filename,
    column_types=None,
    exclude_columns=None,
    include_unaligned=True,
    block_size=None,
    bam_index=None,
):
    if column_types is None:
        column_types = {}
    if exclude_columns is None:
        exclude_columns = []
    exclude_columns = set(exclude_columns)
    if isinstance(bam_filename, pysam.AlignmentFile):
        bam = bam_filename
    else:
        bam = pysam.AlignmentFile(bam_filename, check_sq=False)
    if bam_index is None:
        bam_index = pysam.IndexedReads(bam)
        bam.reset()  # must reset or else bam_index won't include all reads
        bam_index.build()
    aligned_read_names = set()
    bam_types = {}
    if "read_seq" not in exclude_columns:
        bam_types["read_seq"] = pa.string()
    if "read_phred" not in exclude_columns:
        bam_types["read_phred"] = pa.list_(pa.uint8())
    bam_types = {**bam_types, **column_types}
    unaligned_batch_size = None
    new_batch = None
    for batch in iter_gaf(gaf_filename, block_size=block_size):
        if unaligned_batch_size is None:
            # use the size of the first GAF RecordBatch as the unaligned chunk size
            unaligned_batch_size = len(batch)
        bam_columns = None
        name_col = batch.column("name")
        for idx in range(len(name_col)):
            name = name_col[idx].as_py()
            try:
                read = next(bam_index.find(name))
            except:
                raise ValueError(f"BAM missing read '{name}', has alignment in GAF")
            aligned_read_names.add(name)
            if bam_columns is None:
                bam_columns = {col_name: [] for col_name in bam_types.keys()}
            appended_columns = set()
            if "read_seq" in bam_columns:
                bam_columns["read_seq"].append(read.query_sequence)
                appended_columns.add("read_seq")
            if "read_phred" in bam_columns:
                bam_columns["read_phred"].append(np.asarray(read.query_qualities))
                appended_columns.add("read_phred")
            for tag, value, type_ in read.get_tags(with_value_type=True):
                if tag in exclude_columns:
                    continue
                if tag in bam_columns:
                    bam_columns[tag].append(value)
                else:
                    bam_columns[tag] = [None] * idx + [value]
                    bam_types[tag] = pyarrow_type_for_bam(type_, value)
                appended_columns.add(tag)
            for tag in bam_columns.keys() - appended_columns:
                bam_columns[tag].append(None)
        columns = dict(zip(batch.column_names, batch.columns))
        for col_name, col_type in bam_types.items():
            columns[col_name] = pa.array(bam_columns[col_name], col_type)
        new_batch = pa.RecordBatch.from_pydict(columns)
        yield new_batch
    if include_unaligned:
        bam.reset()  # rewind bam file to beginning
        bam_iter = bam.fetch(until_eof=True)
        bam_columns = None
        num_reads = 0
        bam_types = {"name": pa.string(), **bam_types}
        eof = False
        while True:
            for read in bam_iter:
                if read.query_name not in aligned_read_names:
                    if bam_columns is None:
                        bam_columns = {col_name: [] for col_name in bam_types.keys()}
                    appended_columns = set()
                    if "read_seq" in bam_columns:
                        bam_columns["read_seq"].append(read.query_sequence)
                        appended_columns.add("read_seq")
                    if "read_phred" in bam_columns:
                        bam_columns["read_phred"].append(
                            np.asarray(read.query_qualities)
                        )
                        appended_columns.add("read_phred")
                    for tag, value, type_ in read.get_tags(with_value_type=True):
                        if tag in exclude_columns:
                            continue
                        if tag in bam_columns:
                            bam_columns[tag].append(value)
                        else:
                            bam_columns[tag] = [None] * idx + [value]
                            bam_types[tag] = pyarrow_type_for_bam(type_, value)
                        appended_columns.add(tag)
                    for tag in bam_columns.keys() - appended_columns:
                        bam_columns[tag].append(None)
                    num_reads += 1
                    if num_reads == unaligned_batch_size:
                        break
            else:
                eof = True
            if bam_columns is None:
                # no unaligned reads
                return
            columns = {}
            if new_batch is None:
                # all reads are unaligned, so we don't have a GAF schema
                unaligned_types = bam_types.items()
            else:
                unaligned_types = zip(new_batch.schema.names, new_batch.schema.types)
            for col_name, col_type in unaligned_types:
                if col_name in bam_columns:
                    columns[col_name] = pa.array(bam_columns[col_name], col_type)
                else:
                    columns[col_name] = pa.nulls(num_reads, col_type)
            new_batch = pa.RecordBatch.from_pydict(columns)
            yield new_batch
            if eof:
                return
            # start new batch
            bam_columns = None
            num_reads = 0


def read_bam_and_gaf(bam_filename, gaf_filename, **kwargs):
    return pa.Table.from_batches(iter_bam_and_gaf(bam_filename, gaf_filename, **kwargs))


def format_fastx(seqs, phreds=None, names=None):
    if phreds is None:
        format = "fasta"
        phreds = it.repeat(None)
    else:
        format = "fastq"
    if names is None:
        names = it.repeat(None)
    for idx, (name, seq, phred) in enumerate(zip(names, seqs, phreds)):
        if name is None:
            name = f"seq_{idx}"
        if phred is None:
            letter_annotations = None
        else:
            letter_annotations = dict(phred_quality=phred)
        record = SeqRecord(
            Seq(seq),
            id=name,
            description="",
            letter_annotations=letter_annotations,
        )
        yield record.format(format)


def write_fastx(filename, seqs, phreds=None, names=None):
    with open(filename, "w") as f:
        for s in format_fastx(seqs, phreds=phreds, names=names):
            f.write(s)
