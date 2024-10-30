# Overview
This repo provides three things:
1. An end-to-end Nextflow pipeline that takes un-basecalled POD5 (Oxford Nanopore) or basecalled long reads (Oxford Nanopore or PacBio) in BAM or FASTQ format along with a GFA reference and outputs structured tables of GFA-aligned consensus sequences in Arrow or Parquet format for convenient downstream analysis.
2. Python scripts that execute each step of the pipeline. These are designed to work together (i.e., the output of one script can be read in as input to the next script). If a user has a workflow not supported by the Nextflow pipeline (and the user doesn't wish to modify the Nextflow pipeline), the user could save time by using these scripts as building blocks when constructing their custom pipeline.
    - `chunk_seqs.py`: Converts a set of SAM, BAM, or FASTQ input files to a set of SAM, BAM, or FASTQ output files with a target number of output files or target output file size (in bytes or number of sequences).
    - `filter_gfa.py`: Outputs a GFA file containing a subgraph of the graph specified by the input GFA. You can specify which segments to include/exclude by name or by name prefix.
    - `join_gaf.py`: Takes a GAF alignment produced by GraphAligner and joins it to the original reads (in BAM, Arrow, or Parquet format), outputting an Arrow or Parquet table.
    - `find_duplex_pairs.py`: Given a barcode GFA, GAF alignment, and BAM-formatted reads, outputs a space-separated text file of read ids that can be used by the Oxford Nanopore `dorado duplex` basecaller. Identifies pairs of reads that share the same barcode alignment that were sequenced by the same pore within a specified maximum time delta.
    - `prepare_reads.py`: Takes a table of aligned reads (Arrow/Parquet) as input and outputs a table with additional columns.
    - `consensus.py`: If `--skip-consensus` and `--group X/Y` are given, outputs a table of the Xth read group (of Y total groups) where reads are grouped by their barcode alignment. If `--skip-consensus` is not given, outputs a table of consensus sequences for each read group found in the input reads. Options allow filtering which reads are used to build consensus sequences.
    - `realign.py`: Given a table of consensus sequences aligned to the variants GFA, this uses either the `parasail` (default) or `pywfa` pairwise aligner to realign consensus sequences to the linear reference constructed using the GAF alignment. This pairwise realignment is optional but compared with the GraphAligner output offers better guarantees of alignment optimality and more flexible handling of degenerate bases and nonstandard alignment penalties.
    - `extract_segments.py`: Given a table of consensus sequences aligned to the variants GFA (and optionally realigned), this outputs a struct column. For each segment (by default all segments in the GFA, options allow including/excluding segments) the output table contains the columns specifying the number of matches/mismatches/insertions/deletions across that segment, the slice of the CIGAR string corresponding to each segment, the slice of the consensus sequence corresponding to each segment, and the variant name corresponding to that segment (if applicable).
3. Python libraries that provide the core functionality exposed by the scripts. Of particular note:
    - `io.iter_bam_and_gaf`: Joins reads in BAM format with a GAF alignment.
    - `io.iter_gaf`: Parses GAF alignments (tested with GraphAligner's GAF output).
    - `processing.find_duplex_pairs`: Implements a polars query to identify Oxford Nanopore duplex read pairs with matching barcode alignments that passed through the same pore within a specified maximum time delta.
    - `align.pairwise_align`: Generalized interface to pairwise aligners `parasail` and `pywfa`. Makes handling degenerate bases with `parasail` easy, returns a uniform CIGAR representation.
    - `processing.pairwise_align_df_to_path`: Convenience function for running `align.pairwise_align` on every row in a polars dataframe.
    - `consensus.poa`: Generalized interface to consensus generation using `abpoa` or `spoa`.
    - `cigar.cut_cigar`: The algorithm underlying the `extract_segments.py` script. Partitions consensus sequences and CIGAR strings into subcolumns for each GFA segment. Outputs a format that makes downstream analysis far easier. See description of `EXTRACT_SEGMENTS` step [below](#extract-segments) for details.
    - `processing.cut_cigar_df`: Convenience function for running `cigar.cut_cigar` on every row in a polars dataframe.

# GFAs
[GFAs](https://gfa-spec.github.io/GFA-spec/GFA1.html) are graph reference genomes that specify the expected sequence content. (Not to be confused with [GAFs](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf), which are the SAM-style alignment files produced by GraphAligner given FASTQ reads and a GFA reference as inputs.) The pipeline requires two GFAs as input: `grouping_gfa` and `variants_gfa`. The grouping GFA contains segments used to group reads; these groups each result in a consensus sequence. These consensus sequences are then aligned against the variants GFA. Typically, the grouping GFA contains a barcode and the variants GFA contains both the barcode and the construct of interest. Note that all scripts take options that allow subsetting GFA segments by name or by prefix (`--include`, `--exclude`, `--include-prefix`, `--exclude-prefix`). In principle one could use the same GFA file for both grouping and variants, and just add options to exclude all the construct segments during the grouping step, but it is generally easier to use two GFA files (where the grouping GFA file is a subgraph of the variants GFA). This can be done manually; alternatively, the `filter_gfa.py` script may be used to write a filtered GFA file applying the given include/exclude options. Using the same GFA for both `grouping_gfa` and `variants_gfa` may be useful in cases where you wish to group reads using both the barcode and variant segments (e.g., if you want to detect barcode collisions).

The pipeline assumes the convention that different sequence variants of the same logical sequence region may be named `my_segment=variant1`, `my_segment=variant2`, etc. These segments should be mutually exclusive (the topology of the GFA should be such that at most one may appear in any sequence, this is not checked). In the output of `extract_segments.py`, the columns corresponding to these segments will be collapsed into a single set of columns (i.e., the number of mismatches for that sequence region will be recorded in `my_segment|mismatches` instead of `my_segment=variant1|mismatches` or `my_segment=variant2|mismatches`). The additional column `my_segment|variant` will specify which variant appeared (`variant1`, `variant2`, etc). If this behavior is not desired, you can pass `--no-variant-sep` to `extract_segments.py`; it also accepts a `--variant-sep =` option to specify an alternative separator string.

Because `extract_segments.py` slices sequence regions for each GFA segment, it is convenient to have each logically-meaningful variable sequence region (e.g., promoter, RBS) specified as its own GFA segment. Note that GraphAligner seeding will only work on longer segments (e.g., CDSes), so you should avoid GFAs that are only made up of lots of extremely short segments with no long segments, but this should not be a practical limitation.

# Pipeline
1. Different paths are taken depending on what type of input is given:
    - FASTQ input (`basecall: false` and `fastq_input` is non-null)
        1. `SAMTOOLS_IMPORT`: Convert input FASTQ to BAM.
        1. If `publish_bam` is `true`, symlink BAM to `path/to/data/output/RUN_PATH/bam`.
        1. `GRAPHALIGNER_GROUPING`: Align reads to grouping GFA.
    - BAM input (`basecall: false` and `bam_input` is non-null)
        1. `SAMTOOLS_FASTQ`: Convert input BAM to FASTQ.
        1. If `publish_fastq` is `true`, symlink FASTQ to `path/to/data/output/RUN_PATH/fastq`.
        1. `GRAPHALIGNER_GROUPING`: Align reads to grouping GFA.
    - Raw Oxford Nanopore reads as POD5 input (`basecall: true` and `pod5_input` is non-null)
        1. `POD5_MERGE`: Optionally combine input POD5 files into a smaller number of larger POD5 files.
        1. `POD5_VIEW_AND_SUBSET`: Reorganizes raw reads into multiple output files grouped by read properties. By default, this outputs a POD5 file for each nanopore channel. To parallelize this I/O- and memory-intensive operation, a `POD5_VIEW_AND_SUBSET` task is run for each input POD5 file (or output of `POD5_MERGE`).
        1. `POD5_MERGE2`: Multiple `POD5_VIEW_AND_SUBSET` tasks will result in multiple POD5 files for each channel, this merges those into a single POD5 file per channel.
        1. Different basecalling methods are available:
            - Simplex only (`duplex: false`)
                1. `DORADO_DOWNLOAD`: Download simplex dorado basecalling model.
                1. `DORADO_BASECALLER`: Simplex basecall.
                1. If `publish_bam` is `true`, symlink BAM to `path/to/data/output/RUN_PATH/bam`.
            - Duplex with dorado pairing (`duplex: true, use_dorado_duplex_pairing: true`)
                1. `DORADO_DOWNLOAD`: Download simplex dorado basecalling model.
                1. `DORADO_DOWNLOAD2`: Download duplex dorado basecalling model.
                1. `DORADO_DUPLEX`: Duplex basecall.
                1. If `publish_bam` is `true`, symlink BAM to `path/to/data/output/RUN_PATH/bam`.
                1. If `publish_fastq` is `true`, symlink FASTQ to `path/to/data/output/RUN_PATH/fastq`.
                1. `GRAPHALIGNER_VARIANTS`: Align reads to grouping GFA.
            - Duplex with non-dorado pairing (`duplex: true, use_dorado_duplex_pairing: false`)
                1. `DORADO_DOWNLOAD`: Download simplex dorado basecalling model.
                1. `DORADO_DOWNLOAD2`: Download duplex dorado basecalling model.
                1. `DORADO_BASECALLER`: Simplex basecall.
                1. If `publish_bam_simplex` is `true`, symlink BAM to `path/to/data/output/RUN_PATH/bam_simplex`.
                1. `SAMTOOLS_FASTQ`: Convert simplex BAM to FASTQ.
                1. `GRAPHALIGNER_GROUPING`: Align simplex reads to grouping GFA.
                1. `FIND_DUPLEX_PAIRS`: Output a space-separated text file of read ID pairs for valid duplex pairs (a valid duplex pair has two parent reads which align to the same barcode and traversed the same nanopore within a small time delta of each other).
                1. `DORADO_DUPLEX_WITH_PAIRS`: Duplex basecall only those pairs of reads identified by `FIND_DUPLEX_PAIRS`.
                1. `SAMTOOLS_MERGE`: Combine simplex and duplex BAMs.
                1. If `publish_bam` is `true`, symlink BAM to `path/to/data/output/RUN_PATH/fastq`.
                1. `SAMTOOLS_FASTQ_DUPLEX`: Convert duplex BAM to FASTQ.
                1. `CAT_DUPLEX_FASTQ`: Combine simplex and duplex FASTQs.
                1. If `publish_fastq` is `true`, symlink FASTQ to `path/to/data/output/RUN_PATH/fastq`.
                1. `GRAPHALIGNER_GROUPING_DUPLEX`: Align duplex reads to grouping GFA.
                1. `CAT_DUPLEX_GAF`: Combine simplex and duplex GAFs.
1. `JOIN_GAF_GROUPING`: Join reads with alignment against grouping GFA.
1. `PREPARE_READS`: Adds the following columns:
    - A `path` column which is a normalized alignment path (a subset of segments are kept according to the given grouping GFA and segment include/exclude options, and the orientation is normalized; the `reverse_complement` column records the read orientation). The original aligment path is kept as the `full_path` column. This path (specifically, the `path_hash` column) is used for read grouping in the consensus step.
    - A `path_hash` integer column which is a hash of the `path` column. Note that this hash is not stable across polars versions.
    - An `end_to_end` column indicating if GAF alignment path covers all expected segments (specified via grouping GFA and segment include/exclude options).
    - A `grouping_segments` struct column containing matches/mismatches/insertions/deletions/divergence for each grouping segment.
    - If `--max-divergence` is given to `prepare_reads.py`, a `max_divergence` column is added, and it gives the maximum divergence across all segments.
    - If the `qs` BAM tag was not included in input data (it is output by dorado for Oxford Nanopore data, it is removed from PacBio input if the `rq` tag is present because `qs` does not mean quality in PacBio output), it is recomputed as the mean phred score.
    - For Oxford Nanopore duplex data (where the `dx` BAM tag is present), the `is_valid` column is added. Duplex reads with parents whose barcodes match and simplex reads which are not a parent to a valid duplex read.
1. If `publish_prepare_reads` is `true`, symlink Arrow/Parquet output to `path/to/data/output/RUN_PATH/prepare_reads`.
1. Consensus.
    - If `prepare_consensus` is `true` is set, consensus generation proceeds in two steps. This is the default as it may be more memory efficient to specify different memory limits for the preparation (read grouping) step and the consensus generation step.
        1. `PREPARE_CONSENSUS`: Groups reads by alignment path into a set of output files (number specified by `consensus_jobs`). If `--limit-depth N` is given, `read_seq`/`read_phred`/`grouping_segments` columns will be set to `null` for every read after the Nth read in each group (this is to reduce storage/memory usage when a small number of groups has a disproportionate number of reads assigned to it, usually these groups are not complete barcodes and will be filtered out anyway, but they are retained in this step for diagnostic purposes). Adds the following columns:
            - `grouping_depth`: Number of reads which share that read's alignment path.
            - `grouping_duplex_depth`: Same, but only counts duplex reads. Only generated in input reads contain `dx` SAM tag.
        1. If `publish_prepare_reads` is `true`, symlink Arrow/Parquet output to `path/to/data/output/RUN_PATH/prepare_reads`.
        1. `CONSENSUS_PREPARED`: For each input file, compute consensus sequences for each read group in that file. If `--min-depth` or `--min-duplex-depth`, only computes consensus sequences for groups with at least those depths. If `--limit-depth N` is given, only use the longest `N` sequences to compute consensus sequences, discarding the remainder. Adds the columns `consensus_depth` (and `consensus_duplex_depth` if the `dx` SAM tag appears in input) recording the depth that was used for computing the consensus sequence. If `--max-divergence` is given, only use reads with alignments with at most that maximum divergence across all segments. `--max-length L` filters out all reads longer than `L`.
    - `CONSENSUS`: Do both `PREPARE_CONSENSUS` and `CONSENSUS_PREPARED` in a single step.
1. If `publish_consensus` is `true`, symlink Arrow/Parquet output to `path/to/data/output/RUN_PATH/consensus`.
1. `GRAPHALIGNER_VARIANTS`: Align consensus sequences to variants GFA.
1. `JOIN_GAF_VARIANTS`: Join consensus sequences with output of `GRAPHALIGNER_VARIANTS`. This will output alignment path as `variants_path` and alignment CIGAR string as `cg`.
1. If `publish_join_gaf_variants` is `true`, symlink Arrow/Parquet output to `path/to/data/output/RUN_PATH/join_gaf_variants`.
1. `REALIGN`: Given a table of consensus sequences aligned to the variants GFA, `realign.py` uses either the `parasail` (default) or `pywfa` pairwise aligner to realign consensus sequences to the linear reference constructed using the GAF alignment. This pairwise realignment is optional but compared with the GraphAligner output offers better guarantees of alignment optimality and more flexible handling of degenerate bases and nonstandard alignment penalties. If `realign` is `true`, realignment score and CIGAR string will be stored in `realign_score`, `realign_cg` columns, respectively.
1. If `publish_realign` is `true`, symlink Arrow/Parquet output to `path/to/data/output/RUN_PATH/realign`.
1. <a id="extract-segments"></a>`EXTRACT_SEGMENTS`: Uses the consensus sequence, alignment path, and alignment CIGAR string to define a large number of columns that are convenient for downstream analysis. If `realign` is `false`, `--cigar-col cg` must be included in `extract_segments_args` (it defaults to `realign_cg`). Adds a `is_primary_alignment` column indicating if that row corresponds to the primary alignment of a given consensus sequence to the variants GFA (non-primary alignments should be filtered out for most downstream analysis). It also adds the following columns for each segment `seg`:
    - `seg|seq`: The slice of the consensus sequence that aligns to segment `seg`.
    - `seg|cigar`: The slice of the CIGAR string that corresponds to the alignment across segment `seg`.
    - `seg|variant` (if applicable): If there are mutually exclusive segments `seg=variant1`, `seg=variant2`, and so forth, this specifies which of those variants (`variant1`, `variant2`, etc.) appeared in the alignment.
    - `seg|matches`: The number of matches in the CIGAR string for this segment.
    - `seg|mismatches`: The number of mismatches in the CIGAR string for this segment.
    - `seg|insertions`: The number of insertions in the CIGAR string for this segment.
    - `seg|deletions`: The number of deletions in the CIGAR string for this segment.
1. If `publish_extract_segments` is `true`, symlink Arrow/Parquet output to `path/to/data/output/RUN_PATH/extract_segments`.

# Installation
Conda environment creation takes a while, so do this in an interactive job on HMS O2.
```bash
git clone git@github.com:paulssonlab/barcoded_amplicon_pipeline.git
cd barcoded_amplicon_pipeline
mamba create -f environment.yml -n barcoded_amplicon_pipeline
mamba activate barcoded_amplicon_pipeline
pth_file="$(python -c 'import site; print(site.getsitepackages()[0])')/barcoded_amplicon_pipeline.pth"
echo "$PWD/paulssonlab/src" > "$pth_file"
```

This creates a conda environment called `barcoded_amplicon_pipeline` (you can call it whatever you want by changing the argument to `-n`). Whenever you want to use the pipeline, you will first need to activate it with `mamba activate barcoded_amplicon_pipeline`. In the below instructions, `$PIPELINE_ROOT` refers to the directory into which this repo was cloned.

# Usage
1. Copy the input data (FASTQ, BAM, or POD5) to a scratch directory (on HMS O2, it should be on the scratch filesystem).
2. FASTQ and BAM output from Oxford Nanopore or PacBio instruments need to be rechunked (into more, smaller files). In an interactive job, with the `barcoded_amplicon_pipeline` conda environment activated, run `python $PIPELINE_ROOT/sequencing/bin/chunk_seqs.py --size 2147483648 "~/path/to/data/fastq_pass/*.fastq.gz" ~/path/to/data/fastq_chunked`. 2147483648 bytes corresponds to 2 GiB, which seems to work well. This may take 30-60 min for a PromethION-sized sequencing run.
3. Create a [GFA](https://gfa-spec.github.io/GFA-spec/GFA1.html). See examples in `examples/gfas`.
4. Make a copy of the configuration run directory `examples/run`.
5. In that run directory, edit `nextflow.config` so that `params.root` points to the scratch directory containing the data and `samples.toml` to specify the desired input and output (see [Configuration](#configuration)). If you rechunked FASTQ or BAM as described above, you should be sure to specify the rechunked (not original) FASTQ or BAM as input to the pipeline.
6. `cd` into the run directory and start the pipeline with the command (modifying `PIPELINE_ROOT` to point to where the pipeline repo is cloned): `PIPELINE_ROOT=/home/$USER/barcoded_amplicon_pipeline MAMBA_ALWAYS_YES=true time nextflow run main.nf -profile o2 -with-report -with-timeline > out.log`
7. Open another terminal (or tmux) window and periodically run `less out.log` (scroll to the bottom with shift-G) to monitor progress. Additionally, it may be helpful to monitor SLURM jobs with `squeue --me|less`. If the pipeline exits prematurely, see [troubleshooting](#troubleshooting).
8. Optionally, you may use the `notebooks/Eaton.ipynb` notebook to convert the output of the `extract_segments` step (the final step) into a Daniel Eaton format codebook.

Pipeline output for each step is written to the directory `/path/to/data/output/RUN_PATH/STEP`. By default, `RUN_PATH` is `default`. `STEP` is one of `POD5`, `bam_simplex`, `bam`, `fastq`, `prepare_reads`, `prepare_consensus`, `consensus`, `join_gaf_variants`, `realign`, and `extract_segments`. Note that the output files contained are symlinks into `/path/to/data/work`. You must include the `-L` option when copying: e.g., archive the pipeline output with a command like `cp -RL /path/to/data/output/default/extract_segments /path/to/archival/storage/my_run_extract_segments`. Note that if you `rm -rf /path/to/data/work` to clear pipeline intermediate data, this will **delete pipeline output** as well. It is recommended that the desired output is copied to long-term storage upon successful completion of the pipeline. If you run nextflow with the `-report` and `-timeline` options, nextflow will output files with names like `report-20241029-36282066.html` and `timeline-20241029-36282066.html` into the run directory. One contains a timeline showing the execution schedule of each task. The other contains quantile plots of memory usage/execution time for each step. It also contains an interactively-filterable table of tasks, statistics about each task, and each task's working directory (under `path/to/data/work`). This last is is helpful for debugging, as you may want to `cd` into one of those task working directories and re-run a computation with `bash .command.sh` to validate a code change before rerunning the pipeline.

# Configuration

When setting up a new run, make a copy of the `examples/run` directory. This directory should contain three files:
- `main.nf`: Boilerplate, does not need to be modified.
- `nextflow.config`: Edit so that `params.root` points to a directory containing input data. By default, Nextflow will store its intermediate files under `${params.root}/work` and pipeline output will be symlinked under `${params.root}/output/RUN_PATH/STEP_NAME`.
- `samples.toml`: Specifies the operation of the pipeline. One of the input parameters (e.g., `pod5_input`, `fastq_input`, `bam_input`) should be set to the correct glob expression.

## `samples.toml`
One or more pipeline runs may be defined in a `samples.toml` file. The `samples.toml` mini-language uses thet [TOML](https://toml.io/) syntax and allows for flexible processing of multiple datasets (“samples”) with multiple different parameter sets each.

A `samples.toml` file can contain:
- Zero or more `[[samples]]` blocks to specify multiple datasets (“samples”) to be processed.
- Zero or more `[[params]]` blocks to specify multiple parameter sets (“param set”) that are to be applied to each sample.
- An optional `[defaults]` block which specifies default parameters (can be overridden by samples and param sets).
- An optional top-level string named `tsv` which contains a table of tab-separated values. The rows of the tsv table become the list of samples and the columns the different parameter values (the first line of the tsv is interpreted as a column header).

The pipeline is run for each combination of sample and param set (i.e., the set of pipeline runs is the Cartesian product of param sets and samples). The parameters given in a `[[samples]]` block override those given in a `[[params]]` blocks, which in turn overrides those given in the `[defaults]` block. Samples are expected to have a unique `name` parameter, and param sets are expected to have a unique `run_path` parameter. In that case, the output directory of a given run will be `${params.root}/output/${run_path}/${name}/STEP_NAME`. If no samples are defined, the output directory will be `${params.root}/output/${run_path}/STEP_NAME`. If no samples nor param sets are defined (only the `[defaults]` block), the output directory will be `${params.root}/output/defaults/STEP_NAME`. String parameter values can use [Groovy-style string interpolation](https://groovy-lang.org/syntax.html#_string_interpolation) to automatically make substitutions using other parameter values.

For examples of advanced usage, see the `examples/samples.toml` directory.

The following simple examples illustrate typical usage.

An example for processing PromethION FASTQ output:
```toml
[defaults]
fastq_input = "fastq_chunked/*.fastq.gz"
basecall = false
gfa_grouping = "/home/jqs1/scratch/sequencing/sequencing_references/barcode.gfa"
gfa_variants = "/home/jqs1/scratch/sequencing/sequencing_references/pLIB492-501.gfa"
consensus_args = "--method abpoa --no-phred-output --max-length 10000 --max-divergence 0.3"
consensus_jobs = 400
consensus_jobs_per_align_job = 8
find_duplex_pairs_args = "-x UNS9,BC:T7_prom,BC:UMI:upstream,BC:UMI,BC:UMI:downstream,BC:spacer2,BC:term:T7,BC:term:T7hyb10,JUNC10_UNS10"
prepare_reads_args = "-x UNS9,BC:T7_prom,BC:UMI:upstream,BC:UMI,BC:UMI:downstream,BC:spacer2,BC:term:T7,BC:term:T7hyb10,JUNC10_UNS10"
```

Or for basecalling PromethION POD5:
```toml
[defaults]
pod5_input = "pod5_pass/*.pod5"
pod5_chunk = true
pod5_split = true
basecall = true
duplex = true
use_dorado_duplex_pairing = false
dorado_model = "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
dorado_duplex_model = "dna_r10.4.1_e8.2_5khz_stereo@v1.3"
gfa_grouping = "/home/jqs1/scratch/sequencing/sequencing_references/barcode.gfa"
gfa_variants = "/home/jqs1/scratch/sequencing/sequencing_references/pLIB492-501.gfa"
consensus_args = "--method abpoa --no-phred-output --max-length 10000 --max-divergence 0.3"
consensus_jobs = 400
consensus_jobs_per_align_job = 8
find_duplex_pairs_args = "-x UNS9,BC:T7_prom,BC:UMI:upstream,BC:UMI,BC:UMI:downstream,BC:spacer2,BC:term:T7,BC:term:T7hyb10,JUNC10_UNS10"
prepare_reads_args = "-x UNS9,BC:T7_prom,BC:UMI:upstream,BC:UMI,BC:UMI:downstream,BC:spacer2,BC:term:T7,BC:term:T7hyb10,JUNC10_UNS10"
```

These examples include a comma-delimited list of segments specified with the `-x` (exclude) option to `find_duplex_pairs.py` and `prepare_reads_args.py`. This excludes the constant (non-variable) segments in the barcode, making it such that grouping is only done using the `BC:BIT0=0`,...,`BC:BIT29=1` segments. Without this, one read that has the same barcode as another read except it truncates right after the barcode so does not align to one of the constant segments flanking the barcode bits would be placed in a different read group. This is likely not desirable, hence these exclusions. Note that `consensus.py --skip-consensus` (run during the `PREPARE_CONSENSUS` step) actually performs the read grouping, but `prepare_reads.py` (during the `PREPARE_READS` step) performs the orientation normalization and filtering of the alignment path which determines how reads will be grouped. In general, you likely want to supply the same segment include/exclude options to both `find_duplex_pairs.py` and `prepare_reads_args.py`, though `find_duplex_pairs_args` is only used for runs with the parameters `basecall: true, duplex: true, use_dorado_duplex_pairing: false`.

## Parameters
The pipeline can be configured with the following parameters:
<dl>
<dt><code>pod5_input</code><dt>
<dd>Path of POD5 files relative to the <code>params.root</code> directory. Can contain wildcards (i.e., <code>"pod5_pass/*.pod5"</code>).</dd>
<dt><code>bam_input</code><dt>
<dd>Path of BAM files relative to the <code>params.root</code> directory. Can contain wildcards (i.e., <code>"bam/*.bam"</code>).</dd>
<dt><code>bam_simplex_for_duplex_input</code><dt>
<dd>Path of BAM files relative to the <code>params.root</code> directory. Can contain wildcards (i.e., <code>"bam_simplex/*.bam"</code>). If <code>bam_simplex_for_duplex_input</code> and <code>pod5_input</code> are both given and <code>duplex: true, use_dorado_duplex_pairing: false</code> is specified, duplex basecalling is performed without repeating simplex basecalling.</dd>
<dt><code>fastq_input</code><dt>
<dd>Path of FASTQ files relative to the <code>params.root</code> directory. Can contain wildcards (i.e., <code>"fastq_pass/*.fastq"</code>).</dd>
<dt><code>prepare_reads_input</code><dt>
<dd>Path of Arrow/Parquet files (output from the <code>PREPARE_READS</code> step) relative to the <code>params.root</code> directory. Can contain wildcards (i.e., <code>"prepare_reads/*.arrow"</code>).</dd>
<dt><code>prepare_consensus_input</code><dt>
<dd>Path of Arrow/Parquet files (output from the <code>PREPARE_CONSENSUS</code> step) relative to the <code>params.root</code> directory. Can contain wildcards (i.e., <code>"prepare_consensus/*.arrow"</code>).</dd>
<dt><code>consensus_tabular_input</code><dt>
<dd>Path of Arrow/Parquet files (output from the <code>CONSENSUS</code> or <code>CONSENSUS_PREPARED</code> steps) relative to the <code>params.root</code> directory. Can contain wildcards (i.e., <code>"consensus/*.arrow"</code>). Both this and <code>consensus_fasta_input</code> must be given.</dd>
<dt><code>consensus_fasta_input</code><dt>
<dd>Path of FASTA files (output from the <code>CONSENSUS</code> or <code>CONSENSUS_PREPARED</code> steps) relative to the <code>params.root</code> directory. Can contain wildcards (i.e., <code>"consensus/*.fasta"</code>). Both this and <code>consensus_fasta_input</code> must be given.</dd>
<dt><code>realign_input</code><dt>
<dd>Path of Arrow/Parquet files (output from the <code>REALIGN</code> step) relative to the <code>params.root</code> directory. Can contain wildcards (i.e., <code>"realign/*.arrow"</code>).</dd>
<dt><code>basecall</code> (default: <code>true</code>)<dt>
<dd>Whether to basecall Oxford Nanopore POD5 input. If <code>true</code>, <code>POD5_input</code> must be specified.</dd>
<dt><code>duplex</code> (default: <code>true</code>)<dt>
<dd>Whether to use Oxford Nanopore duplex basecalling.</dd>
<dt><code>use_dorado_duplex_pairing</code> (default: <code>false</code>)<dt>
<dd>If <code>true</code>, <code>dorado duplex</code> is used to duplex basecall POD5 input. The majority of duplex pairs will be invalid. <code>prepare_reads.py</code> is used to flag which duplex reads are valid (i.e., both parent reads align to the same barcode); invalid duplex pairs will be filtered out when generating consensus sequences. If <code>false</code>, <code>dorado basecaller</code> is used to simplex basecall POD5 input. These simplex reads are aligned, and <code>find_duplex_pairs.py</code> is used to identify valid candidate duplex pairs. Then <code>dorado duplex</code> is used to duplex basecall only those pairs. Duplex and simplex reads and their alignments are concatenated for subsequent processing. Note that </dd>
<dt><code>dorado_model</code><dt>
<dd>The name of the dorado basecalling model to use for simplex basecalling (<code>dorado model --list</code> prints a list of available models). Must be specified if simplex basecalling is to be performed.
<dt><code>dorado_duplex_model</code><dt>
<dd>The name of the dorado basecalling model to use for duplex basecalling (<code>dorado model --list</code> prints a list of available models). Must be specified if duplex basecalling is to be performed.
<dt><code>pod5_chunk</code> (default: <code>false</code>)<dt>
<dd>Whether POD5 chunking happens before and after POD5 splitting (<code>POD5_VIEW_AND_SUBSET</code>).</dd>
<dt><code>pod5_chunk_before_split</code>, <code>pod5_chunk_after_split</code> (default: none)<dt>
<dd>Can be used to individually set whether POD5 chunking happens before and/or after POD5 splitting (<code>POD5_VIEW_AND_SUBSET</code>). By default, these are both set to the value of <code>pod5_chunk</code>.</dd>
<dt><code>pod5_chunk_bytes</code> (default: <code>10737418240</code>, i.e., 10 GiB)<dt>
<dd>If any of <code>pod5_chunk</code>, <code>pod5_chunk_before_split</code>, or <code>pod5_chunk_after_split</code> are <code>true</code>, specifies a target POD5 file size. Input POD5 files will be grouped before combining such that the resulting POD5 files are each at least this size. Note that if the input POD5 files are larger than this size, the pipeline will not split input files to reach this target size.</dd>
<dt><code>pod5_chunk_files</code><dt>
<dd>If any of <code>pod5_chunk</code>, <code>pod5_chunk_before_split</code>, or <code>pod5_chunk_after_split</code> are <code>true</code>, specifies how many input POD5 files will be combined together into one file. If this is set, <code>pod5_chunk_bytes</code> must be set to <code>null</code>.</dd>
<dt><code>dorado_job_bytes</code> (default: <code>10737418240</code>, i.e., 10 GiB)<dt>
<dd>Group POD5 files so that each basecalling job takes at least 10 GiB of input POD5 data.</dd>
<dt><code>dorado_jobs</code><dt>
<dd>This sets the number of parallel dorado jobs. If this is set, <code>dorado_job_bytes</code> must be set to <code>null</code>.</dd>
<dt><code>pod5_split</code> (default: <code>false</code>)<dt>
<dd>Whether to split POD5 files into groups of files where groups are determined by read properties. Oxford Nanopore <a href="https://github.com/nanoporetech/dorado?tab=readme-ov-file#improving-the-speed-of-duplex-basecalling">recommends</a> splitting by nanopore channel to improve duplex basecalling performance on slow filesystems.</dd>
<dt><code>pod5_split_by</code> (default: <code>["channel"]</code>)<dt>
<dd>Which properties of the POD5 raw reads are used to group reads. Setting this to <code>["channel", "mux"]</code> allows splitting by both channel and mux if a higher degree of parallelization is desired (this is unlikely to be necessary).</dd>
<dt><code>realign</code> (default: <code>true</code>)<dt>
<dd>Whether to run the <code>REALIGN</code> step (see <a href="#pipeline">above</a>).</dd>
<dt><code>publish_pod5</code>, <code>publish_bam_simplex</code>, <code>publish_bam</code>, <code>publish_fastq</code>, <code>publish_prepare_reads</code>, <code>publish_prepare_consensus</code>, <code>publish_consensus</code>, <code>publish_join_gaf_variants</code>, <code>publish_realign</code>, <code>publish_extract_segments</code> (default: <code>true</code>)<dt>
<dd>Whether to symlink the output of each step to `path/to/data/output/RUN_PATH/STEP_NAME`.</dd>
<dt><code>output</code> (default: <code>"extract_segments"</code>)<dt>
<dd>Possible values: <code>"pod5"</code> (output pod5 after chunking/splitting), <code>"basecaller"</code> (output basecalled BAM and FASTQ), <code>"consensus"</code> (output consensus sequences), <code>"extract_segments"</code> (run pipeline to completion).</dd>
<dt><code>samtools_fastq_args</code> (default: <code>"-T \"*\""</code>)<dt>
<dd>Arguments to <code>samtools fastq</code> when converting BAM to FASTQ. The default option preserves all BAM tags.</dd>
<dt><code>samtools_import_args</code> (default: <code>"-T \"*\""</code>)<dt>
<dd>Arguments to <code>samtools import</code> when converting FASTQ to BAM. The default option preserves all BAM tags included in FASTQ sequence id line.</dd>
<dt><code>samtools_merge_args</code> (default: <code>"-c"</code>)<dt>
<dd>Arguments to <code>samtools merge</code> when merging BAM files. The default option combines <code>@RG</code> headers.</dd>
<dt><code>tabular_format</code> (default: <code>"arrow"</code>)<dt>
<dd>The format used for intermediate and final tabular output files, either <code>"arrow"</code> or <code>"parquet"</code>.</dd>
<dt><code>graphaligner_args</code> (default: <code>"-x vg"</code>)<dt>
<dd>Arguments to <code>GraphAligner</code>, see its <a href="https://github.com/maickrau/GraphAligner">documentation</a>.</dd>
<dt><code>find_duplex_pairs_args</code> (default: <code>""</code>)<dt>
<dd>Arguments to <code>find_duplex_pairs.py</code>.</dd>
<dt><code>prepare_reads_args</code> (default: <code>""</code>)<dt>
<dd>Arguments to <code>prepare_reads.py</code>.</dd>
<dt><code>prepare_consensus</code> (default: <code>true</code>)<dt>
<dd>If <code>true</code>, <code>PREPARE_CONSENSUS</code> is used to write read groups into files (all reads from a given group are guaranteed to be written to the same output file). For each of these output files, a <code>CONSENSUS_PREPARED</code> task computes consensus sequences for each read group within that file. If <code>false</code>, the <code>CONSENSUS</code> step groups and compute consensus sequences in a single step. Grouping and computing consensus sequences as separate steps may enable more efficient usage of cluster memory by setting different resource requirements for the two steps.</dd>
<dt><code>consensus_args</code> (default: <code>"--method abpoa --no-phred-output --max-length 10000"</code>)<dt>
<dd>It is highly recommended to set a sensibile <code>--max-length</code> (default: 20000), as Oxford Nanopore data often contains a tail of unexpectedly long reads longer than 100kb (often these are failures of dorado's read-splitting algorithm), and these long sequences cause abpoa/spoa memory usage to spike, often crashing the pipeline. <code>abpoa</code> is faster but may be buggier than <code>spoa</code>.</dd>
<dt><code>prepare_consensus_args</code><dt>
<dd>For the same reason just described, it is highly recommended to set a sensible <code>--max-length</code> (default: 20000), as not doing so can significantly impact the execution time of <code>PREPARE_CONSENSUS</code>, to the point of having those tasks hit their time limit.</dd>
<dt><code>consensus_jobs</code> (default: <code>400</code>)<dt>
<dd>The number of jobs that are used for computing consensus sequences (equivalently, the number of output files for the <code>PREPARE_CONSENSUS</code> step. 400 is recommended for a typical PromethION run, 100 is recommended for a typical MinION run.</dd>
<dt><code>consensus_jobs_per_align_job</code> (default: <code>8</code>)<dt>
<dd>8 is recommended for a typical PromethION run with ~200k barcodes.</dd>
<dt><code>join_gaf_variants_args</code> (default: <code>"--rename-col path grouping_path --rename-col path_hash grouping_path_hash --rename-gaf-col path variants_path"</code>)<dt>
<dd>Arguments to <code>join_gaf.py</code> used when joining consensus sequences with their alignment to the variants GFA. The default arguments keep the <code>path_hash</code> column used for read grouping available, renamed to <code>grouping_path_hash</code>, and makes the alignment path corresponding to the alignment with the variants GFA accessible as <code>variants_path</code>.</dd>
<dt><code>extract_segments_args</code><dt>
<dd>Arguments to <code>extract_segments.py</code>. If <code>realign</code> is <code>false</code>, you must add the option <code>--cigar-col cg</code>.</dd>
</dl>

# Troubleshooting

## Pipeline crashes
If nextflow crashes, look at `out.log` in the run directory and scroll to the bottom. Nextflow should report what task crashed, its working directory, and an exit code. Exit code 137 means a task hit the memory limit, exit code 140 indicates that it hit the time limit. By `cd`-ing into a task's working directory in an interactive job, you can manually try to re-run the task to help identify the issue. If it's a resource limit issue, you can modify the resource limits specified in `paulssonlab/src/paulssonlab/sequencing/modules/scripts.nf` (this shouldn't be necessary). To restart the pipeline after the issue has been addressed, copy the output for the last completed step that preceeded the crash (using `cp -RL`), edit `samples.toml` to specify that copied output as the new input, and then run the pipeline starting from that step. If you want to clear all output and re-run the pipeline from the beginning, run `rm -rf /path/to/data/work` (make sure you copy any outputs from the `path/to/data/output/default/...` directories before you clear the `work` directory, as those outputs live inside the `work` directory). It may be helpful to clean up log files and reports in the run directory by running `rm -rf .nextflow* *.txt *.html *.log`.

For example, let's say that one of the `CONSENSUS_PREPARED` tasks crashed.
```
# we look at the log file and scroll to the end
less out.log
# inside an interactive job,
cd /path/to/data/work/04/40f42f4416de54b6e1f723bf12cd7c
time -v bash .command.sh
# the time -v command prints duration and peak memory usage, which is useful for checking
# that resource limits are appropriate
# let's say this prints out a Python error, we find the bug and fix it
# then we copy the output of the preeceding step so we can start the pipeline there
cp -RL /path/to/data/output/default/prepare_consensus /path/to/data/prepare_consensus
# to prevent the accumulation of large intermediate data files, we clear the work directory
# make sure you copy all outputs you care about first!
rm -rf /path/to/data/work
# now we cd back into our run directory and edit samples.toml
# so instead of fastq_input="fastq_chunked/*.fastq.gz"
# we have prepare_consensus_input="prepare_consensus/*.arrow"
cd /path/to/run
vim samples.toml
# sometimes during debugging you accumulate a lot of log file clutter, so we can clear that
# make sure you don't clear logs you care about!
rm -rf .nextflow* *.txt *.html *.log
# now we start the run again
PIPELINE_ROOT=/home/$USER/barcoded_amplicon_pipeline MAMBA_ALWAYS_YES=true time nextflow run main.nf -profile o2 -with-report -with-timeline > out.log
```

This process is clunkier than it should be, see [below](#nextflow-run--resume).

# Known issues

## Excluding old HMS O2 compute nodes
Running nextflow with the `-profile o2` option will blacklist the compute nodes listed in the text file `$PIPELINE_ROOT/exclude_non_avx2.txt`. This is necessary because the older compute nodes in HMS O2 have old processors that do not support the modern SIMD instruction sets (e.g., AVX2) that are required by the conda-forge builds of polars and pyabpoa. If there's an individual compute node that is causing crashes or you want to blacklist for any reason, you can add its hostname to `$PIPELINE_ROOT/exclude_non_avx2.txt`

## Paulsson Lab GPU node
I have not gotten the pipeline to run successfully on the Paulsson Lab GPU node. It is likely a matter of recompiling polars and/or abpoa (to use the correct SIMD compiler flags) as well as tweaking the dorado resource requirements.

## POD5 chunking
Duplex basecalling is extremely slow when POD5 files are read from high-latency filesystems (i.e., anything except local SSD). Following Oxford Nanopore's [instructions](https://github.com/nanoporetech/dorado), `POD5_VIEW_AND_SUBSET` reorganizes reads from a set of input POD5 files into a set of output POD5 files where each output file corresponds to a Oxford Nanopore channel (it can also split by channel and mux, if desired).

Before early 2024, MinKNOW output a large number of small POD5 files, which made `POD5_VIEW_AND_SUBSET` run slowly. As such, `POD5_MERGE` was used to chunk these files into a slightly smaller number of larger POD5 files (with target size 10 GiB), which significantly improved the performance of `POD5_VIEW_AND_SUBSET`. In early 2024 MinKNOW started outputting much larger POD5 files (average size 30 GiB). The resource requirements of `POD5_VIEW_AND_SUBSET` may need to be increased to cope with these, or else the `POD5_MERGE` step may need to be modified to allow rechunking large POD5 files into a couple smaller POD5 files (using `pod5 subset` and a list of read IDs). Alternatively, `POD5_VIEW_AND_SUBSET` can be skipped entirely if basecalling is done with POD5 files residing on a local fast SSD (I have not tested this because I did not get the pipeline running on the Paulsson Lab GPU node, and the GPU node's SSD was unavailable until an OS upgrade was performed. As of October 2024, I have not been told it has.)

Note that the POD5 chunking/splitting process takes many hours. You may wish to copy the `path/to/data/output/default/pod5` chunked-and-split output and use that as input if you ever wish to rebasecall that dataset in the future (in that case, make sure to set `pod5_chunk` and `pod5_split` to `false` in `samples.toml`).

## Basecalling
[Dorado](https://github.com/nanoporetech/dorado) must be installed if basecalling is to be used. Dorado basecalling models v5.0.0 and newer have significantly longer basecalling runtimes, so resource requirements and chunking settings may need to be adjusted (I have not tested the pipeline with these newer basecalling models). Duplex basecalling has not been tested in a while and may require tweaks to work with the latest version of dorado.

## `nextflow run -resume`
Due to a bug in Nextflow, the `-resume` feature does not work with this pipeline. This means that if some `CONSENSUS_PREPARED` tasks completed but one or more failed and the pipeline aborted, you cannot re-run only those failed tasks using `nextflow run -resume`, you must re-run all `CONSENSUS_PREPARED` tasks. As described [above](#pipeline-crashes), you do this by copying the output for the last completed step that preceeded the crash (using `cp -RL`), edit `samples.toml` to specify that copied output as the new input, and then run the pipeline starting from that step. This means that recovering from pipeline failures can at worst take a couple extra hours compared to what would be the case if `-resume worked`, and thus it is worth submitting a bug report and seeing if Nextflow can address the issue.

Another approach is to replace the use of `map_call_process` with a new custom operator, now that Nextflow supports operator plugins. This could make the implementation simpler, better, and work with the resume machinery.

## Nextflow is pinned to 24.04.4
Use of `@Grab` in `lib/SampleSheetParser.groovy` needs to be replaced, and the code in `lib/` needs to be moved to a Nextflow plugin due to a breaking change in 24.10.0 (see https://github.com/nextflow-io/nextflow/issues/5441).

## Polars memory management
I identified an issue where polars does not release memory (see https://github.com/pola-rs/polars/issues/19497) after collecting a query. My best guess is that it is an interaction between polars, jemalloc, and the extremely old Linux kernel (CentOS 7, Linux 3.10) running on HMS O2. This appears to only a problem during `PREPARE_CONSENSUS` (the most memory-intensive step). I implemented a hacky workaround in `consensus.py` which appears to avoid this issue. There is a comment in `consensus.py` that describes how to remove the workaround if the issue is ever fixed.

## `abpoa` is pinned to a specific 1.5.2 build
See https://github.com/yangao07/abPOA/issues/75#issuecomment-2322807766. It could be that newer builds are only broken on the ancient CPUs of HMS O2; newer builds may work on the Paulsson Lab GPU node or new HMS O2 nodes once they are deployed.
