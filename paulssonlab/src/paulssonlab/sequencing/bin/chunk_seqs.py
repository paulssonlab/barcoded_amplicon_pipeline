#!/usr/bin/env python
import gzip
import itertools as it
import sys
from pathlib import Path

import click
import dnaio
import pysam

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.util import detect_format

FASTQ_FORMATS = ["fastq.gz", "fq.gz", "fastq", "fq"]
FASTA_FORMATS = ["fasta.gz", "fa.gz", "fasta", "fa"]
FORMATS = ["bam", "sam", *FASTQ_FORMATS, *FASTA_FORMATS]


def divide_ceil(n, d):
    return -(n // -d)


def chunk_seqs(
    input_filenames,
    output_dir,
    output_prefix=None,
    input_format=None,
    output_format=None,
    num_seqs=None,
    num_files=None,
    max_size=None,
):
    if sum([num_seqs is not None, num_files is not None, max_size is not None]) != 1:
        raise ValueError(
            "exactly one of num_seqs, num_files, and max_size must be given"
        )
    input_format = detect_format(
        input_format,
        input_filenames[0],
        FORMATS,
    )
    if output_format is None:
        output_format = input_format
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    if output_prefix is not None:
        output_filenames = (
            output_dir / f"{output_prefix}.{idx}.{output_format}" for idx in it.count()
        )
    else:
        output_filenames = (output_dir / f"{idx}.{output_format}" for idx in it.count())
    if input_format in ["sam", "bam"]:
        if output_format == "sam":
            output_mode = "w"
        elif output_format == "bam":
            output_mode = "wb"
        else:
            raise ValueError(f"invalid output format {output_format} for sam/bam input")
        # TODO: dnaio is faster at uBAM reading than pysam/htslib
        # switch if/when dnaio implements write support??
        input_files = (
            pysam.AlignmentFile(filename, check_sq=False)
            for filename in input_filenames
        )
        if num_files is not None:
            input_files = list(input_files)
            total_seqs = sum(
                input_file.count(until_eof=True) for input_file in input_files
            )
            num_seqs = divide_ceil(total_seqs, num_files)
        output_file = None
        header = None
        try:
            for input_file in input_files:
                if header is None:
                    header = input_file.header
                input_file.reset()  # reset file cursor after .count() above
                for seq in input_file:
                    if output_file is None:
                        output_filename = next(output_filenames)
                        output_file = pysam.AlignmentFile(
                            output_filename, output_mode, header=header
                        )
                        seqs_written = 0
                    output_file.write(seq)
                    seqs_written += 1
                    if (num_seqs is not None and seqs_written == num_seqs) or (
                        max_size is not None
                        and output_filename.stat().st_size >= max_size
                    ):
                        output_file.close()
                        output_file = None
        finally:
            if output_file is not None:
                output_file.close()
    elif input_format in [*FASTA_FORMATS, *FASTQ_FORMATS]:
        if output_format not in [*FASTA_FORMATS, *FASTQ_FORMATS]:
            raise ValueError(
                f"invalid output format {output_format} for fasta/fastq input"
            )
        if input_format in FASTA_FORMATS and output_format in FASTQ_FORMATS:
            raise ValueError("cannot output fastq with fasta input")
        input_files = (dnaio.open(filename) for filename in input_filenames)
        if num_files is not None:
            input_files = list(input_files)
            total_seqs = sum(
                1 for filename in input_filenames for _ in dnaio.open(filename)
            )
            num_seqs = divide_ceil(total_seqs, num_files)
        output_file = None
        try:
            for input_file in input_files:
                for seq in input_file:
                    if output_file is None:
                        output_filename = next(output_filenames)
                        output_file = dnaio.open(output_filename, mode="w")
                        seqs_written = 0
                    output_file.write(seq)
                    seqs_written += 1
                    if (num_seqs is not None and seqs_written == num_seqs) or (
                        max_size is not None and output_file._file.tell() >= max_size
                    ):
                        output_file.close()
                        output_file = None
        finally:
            if output_file is not None:
                output_file.close()
    else:
        raise ValueError(f"invalid input format: {input_format}")


@click.command(context_settings={"show_default": True})
@click.option("--input-format", type=click.Choice(FORMATS, case_sensitive=False))
@click.option("--output-format", type=click.Choice(FORMATS, case_sensitive=False))
@click.option("--name", type=str, help="Filename prefix for output files")
@click.option("--seqs", type=int, help="Number of sequences per output file")
@click.option("--files", type=int, help="Number of output files")
@click.option("--size", type=float, help="Max output file size in bytes")
@click.argument("input", type=str, nargs=-1)
@click.argument("output", type=click.Path(file_okay=False))
def cli(
    input,
    output,
    name,
    input_format,
    output_format,
    seqs,
    files,
    size,
):
    if sum([seqs is not None, files is not None, size is not None]) != 1:
        click.error("exactly one of --seqs, --files, --size must be specified")
    chunk_seqs(
        input,
        output,
        name,
        input_format,
        output_format,
        seqs,
        files,
        size,
    )


if __name__ == "__main__":
    cli()
