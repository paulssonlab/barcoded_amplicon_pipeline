#!/usr/bin/env python
import sys
from pathlib import Path

import click
import polars as pl
from gfapy import Gfa

sys.path.append(str(Path(__file__).parents[3]))
from paulssonlab.sequencing.processing import cut_cigar_df
from paulssonlab.sequencing.util import detect_format
from paulssonlab.util.cli import split_delimited_list


def extract_segments(
    gfa_filename,
    input_filename,
    output_filename,
    input_format,
    output_format,
    path_column,
    cigar_column,
    sequence_column,
    phred_column,
    query_start_column,
    query_end_column,
    query_length_column,
    path_start_column,
    path_end_column,
    struct_name,
    keep_full,
    cut_cigar_kwargs,
):
    if not input_filename:
        return  # no op if given no input
    # detect format based on first of potentially many filenames/glob patterns
    input_format = detect_format(
        input_format, input_filename[0], ["arrow", "parquet"], glob=True
    )
    if output_filename:
        output_format = detect_format(
            output_format, output_filename, ["arrow", "parquet"]
        )
    gfa = Gfa.from_file(gfa_filename)
    with pl.StringCache():
        if input_format == "arrow":
            df = pl.concat([pl.scan_ipc(f) for f in input_filename], how="diagonal")
        elif input_format == "parquet":
            df = pl.concat([pl.scan_parquet(f) for f in input_filename], how="diagonal")
        df = cut_cigar_df(
            df,
            gfa,
            path_column=path_column,
            cigar_column=cigar_column,
            sequence_column=sequence_column,
            phred_column=phred_column,
            query_start_column=query_start_column,
            query_end_column=query_end_column,
            query_length_column=query_length_column,
            path_start_column=path_start_column,
            path_end_column=path_end_column,
            struct_name=struct_name,
            keep_full=keep_full,
            cut_cigar_kwargs=cut_cigar_kwargs,
        )
        # TODO: waiting on polars bugs before we can sink_ipc, see realign.py
        df = df.collect()
        if output_format == "arrow":
            df.write_ipc(output_filename)
        elif output_format == "parquet":
            df.write_parquet(output_filename)


@click.command(context_settings={"show_default": True})
@click.option("--gfa", type=click.Path(exists=True, dir_okay=False), required=True)
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
@click.option("--path-col", default="variants_path")
@click.option("--cigar-col", default="realign_cg")
@click.option("--seq-col", default="consensus_seq")
@click.option("--no-seq-col", is_flag=True)
@click.option("--phred-col", default="consensus_phred")
@click.option("--no-phred-col", is_flag=True)
@click.option("--query-start-col", default="query_start")
@click.option("--query-end-col", default="query_end")
@click.option("--query-length-col", default="query_length")
@click.option("--path-start-col", default="path_start")
@click.option("--path-end-col", default="path_end")
@click.option("--pad-alignment", is_flag=True)
@click.option("--segments-struct", default="variants_segments")
@click.option("--no-segments-struct", is_flag=True)
@click.option("--keep-full/--no-keep-full", default=True)
@click.option("-s", "--segments", multiple=True, callback=split_delimited_list)
@click.option("--variant-sep", default="=")
@click.option("--no-variant-sep", is_flag=True)
@click.option("--column-sep", default="|")
@click.option("--indices/--no-indices", default=False)
@click.option("--counts/--no-counts", default=True)
@click.option("--cigars/--no-cigars", default=True)
@click.option("--variants/--no-variants", default=True)
@click.option("--separate-ends/--no-separate-ends", default=True)
@click.argument("input", type=click.Path(exists=True, dir_okay=False), nargs=-1)
@click.argument("output", type=click.Path())
def cli(
    gfa,
    input,
    output,
    input_format,
    output_format,
    path_col,
    cigar_col,
    seq_col,
    no_seq_col,
    phred_col,
    no_phred_col,
    query_start_col,
    query_end_col,
    query_length_col,
    path_start_col,
    path_end_col,
    pad_alignment,
    segments_struct,
    no_segments_struct,
    keep_full,
    segments,
    variant_sep,
    no_variant_sep,
    column_sep,
    indices,
    counts,
    cigars,
    variants,
    separate_ends,
):
    extract_segments(
        gfa,
        input,
        output,
        input_format,
        output_format,
        path_col,
        cigar_col,
        None if no_seq_col else seq_col,
        None if no_phred_col else phred_col,
        None if not pad_alignment else query_start_col,
        None if not pad_alignment else query_end_col,
        None if not pad_alignment else query_length_col,
        None if not pad_alignment else path_start_col,
        None if not pad_alignment else path_end_col,
        None if no_segments_struct else segments_struct,
        keep_full,
        dict(
            # cut_cigar will return empty output if segments is an empty list
            segments=segments or None,
            variant_sep=None if no_variant_sep else variant_sep,
            key_sep=column_sep,
            return_indices=indices,
            return_counts=counts,
            return_cigars=cigars,
            return_variants=variants,
            separate_ends=separate_ends,
        ),
    )


if __name__ == "__main__":
    cli()
