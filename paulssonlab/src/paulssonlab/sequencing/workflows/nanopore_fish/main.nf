import static functions.*

include { POD5_MERGE; POD5_MERGE as POD5_MERGE2; POD5_VIEW_AND_SUBSET } from '../../modules/pod5.nf'
include { DORADO_DOWNLOAD;
          DORADO_DOWNLOAD as DORADO_DOWNLOAD2;
          DORADO_DUPLEX; DORADO_BASECALLER } from '../../modules/dorado.nf'
include { SAMTOOLS_FASTQ } from '../../modules/samtools.nf'
include { GRAPHALIGNER as GRAPHALIGNER_GROUPING; GRAPHALIGNER as GRAPHALIGNER_VARIANTS } from '../../modules/graphaligner.nf'
include { JOIN_GAF as JOIN_GAF_GROUPING;
          JOIN_GAF as JOIN_GAF_VARIANTS;
          PREPARE_READS;
          CONSENSUS;
          REALIGN;
          EXTRACT_SEGMENTS } from '../../modules/scripts.nf'

def GLOBBED_INPUTS = ["bam_input", "fastq_input", "prepare_reads_input", "consensus_tabular_input", "realign_input"]

workflow NANOPORE_FISH {
    take:
    samples_in

    main:
    def DEFAULT_ARGS = [
        tabular_format: "arrow",
        graphaligner_args: "-x dbg",
        prepare_reads_args: "-x UNS9,BC:UPSTREAM,BC:JUNCTION,BC:T7_TERM,BC:SPACER2",
        consensus_args: "--method spoa --no-phred-output --min-depth 3",
        consensus_jobs: 1000,
        consensus_jobs_per_align_job: 10,
        join_gaf_variants_args: "--reads-prefix read_group_ --gaf-prefix consensus_",
        extract_segments_args: "--path-col consensus_path --cigar-col realign_cg",
    ]
    samples_in
    .map {
        if (it.get("basecall") && !it.get("pod5_input")) {
            throw new Exception("pod5_input must be specified to basecall")
        }
        if (it.get("bam_input") && it.get("fastq_input")) {
            throw new Exception("cannot specify both bam_input and fastq_input")
        }
        if (!["basecall", *GLOBBED_INPUTS].collect { k -> it.get(k) }.any()) {
            throw new Exception("one of bam_input or fastq_input required if not basecalling")
        }
        it
    }
    .map {
        it = [*:DEFAULT_ARGS, *:it]
        it = [
            gfa_grouping: it.get("gfa"),
            gfa_variants: it.get("gfa"),
            graphaligner_grouping_args: it.get("graphaligner_args"),
            graphaligner_variants_args: it.get("graphaligner_args"),
            join_gaf_grouping_args: it.get("join_gaf_args"),
            join_gaf_variants_args: it.get("join_gaf_args"),
            *:it
        ]
        if (it.getOrDefault("align", true)) {
            if (!it.get("gfa_grouping")) {
                throw new Exception("gfa or gfa_grouping must be specified if aligning")
            }
            if (!it.get("gfa_variants")) {
                throw new Exception("gfa or gfa_variants must be specified if aligning")
            }
        }
        it
    }
    .branch {
        pod5: it.get("basecall")
        realign: it.get("realign_input")
        consensus: it.get("consensus_tabular_input")
        prepare_reads: it.get("prepare_reads_input")
        bam: it.get("bam_input")
        fastq: it.get("fastq_input")
    }
    .set { ch_input_type }
    ch_input_type.pod5.branch {
        // chunk pod5 files UNLESS:
        // map contains the key pod5_chunk and it is false
        // OR both pod5_chunk_bytes and pod5_chunk_files are falsy
        yes: !(!it.getOrDefault("pod5_chunk", true)
                || (!it.get("pod5_chunk_bytes")
                    && !it.get("pod5_chunk_files")))
        no: true
    }
    .set { ch_do_pod5_chunk }
    ch_do_pod5_chunk.yes.map {
        [*:it, pod5_input_chunked: chunk_files(it["pod5_input"],
                                               it.get("pod5_chunk_files"),
                                               it.get("pod5_chunk_bytes"))
        ]
    }
    .set { ch_pod5_chunked }
    map_call_process(POD5_MERGE,
                     ch_pod5_chunked,
                     ["pod5_input_chunked", "pod5_merge_args"],
                     [id: { pod5, meta -> "${pod5[0].baseName}_merged" }],
                     "pod5_input_chunked",
                     ["pod5_to_split"],
                     "_POD5_MERGE") { pod5, meta -> [pod5] }
    .set { ch_pod5_merged }
    ch_pod5_merged.mix(ch_do_pod5_chunk.no.map {
            [*:it, pod5_to_split: it.get("pod5_input")]
        } )
    .branch {
        yes: it.get("pod5_split")
        no: true
    }
    .set { ch_do_pod5_split }
    ch_do_pod5_split.yes.map {
        def pod5_split_by = it.get("pod5_split_by") ?: ["channel"]
        def pod5_split_chunks = it.get("pod5_split_chunks") ?: 400
        [*:it,
            pod5_split_by: pod5_split_by,
            pod5_view_args: "--include \"${['read_id', *pod5_split_by].join(',')}\" ${it.get('pod5_view_args') ?: ''}",
            pod5_subset_args: "--columns ${pod5_split_by.join(' ')} ${it.get('pod5_subset_args') ?: ''}"]
    }
    .set { ch_to_pod5_split }
    map_call_process(POD5_VIEW_AND_SUBSET,
        ch_to_pod5_split,
        ["pod5_to_split", "pod5_view_args", "pod5_subset_args"],
        [id: { read_list, meta -> read_list.baseName }],
        "pod5_to_split",
        ["pod5_split"],
        "_POD5_VIEW_AND_SUBSET") { pod5_to_split, meta -> [pod5_to_split] }
        .set { ch_pod5_split }
    ch_pod5_split.map { meta ->
        [*:meta, pod5_split_to_chunk: meta.pod5_split.flatten()
                                                        .groupBy { it.baseName }
                                                        .values()
                                                        .collect { it.sort() }
                                                        .sort() { it[0] }
        ]
    }
    .set { ch_to_pod5_split_merge }
    ch_to_pod5_split_merge.branch {
        // chunk pod5 files UNLESS:
        // map contains the key pod5_chunk and it is false
        // OR both pod5_chunk_bytes and pod5_chunk_files are falsy
        yes: !((it.containsKey("pod5_chunk") && !it["pod5_chunk"])
                                    || (!it.get("pod5_chunk_bytes")
                                        && !it.get("pod5_chunk_files")))
        no: true
    }
    .set { ch_do_pod5_split_chunk }
    ch_do_pod5_split_chunk.yes.map {
        [*:it, pod5_split_to_merge: chunk_files(it["pod5_split_to_chunk"],
                                                it.get("pod5_chunk_files"),
                                                it.get("pod5_chunk_bytes"))
        ]
    }
    .set { ch_pod5_split_chunked }
    ch_pod5_split_chunked.mix(ch_do_pod5_split_chunk.no.map {
        [*:it, pod5_split_to_merge: it["pod5_split_to_chunk"]]
    } )
    map_call_process(POD5_MERGE2,
                     ch_pod5_split_chunked,
                     ["pod5_split_to_merge", "pod5_merge_args"],
                     [id: { pod5, meta -> "${pod5[0].baseName}_merged" }],
                     "pod5_split_to_merge",
                     ["pod5_split_merged"],
                     "_POD5_MERGE2") { pod5, meta -> [pod5] }
    .set { ch_pod5_split_merged }
    ch_pod5_split_merged
    .map {
        [*:it, pod5: it.get("pod5_split_merged")]
    }
    .mix(ch_do_pod5_split.no.map { [*:it, pod5: it.get("pod5_to_split")] })
    .set { ch_pod5 }
    ch_pod5.subscribe {
        if (it.get("publish_pod5")) {
            def output_dir = file_in_dir(it.output_run_dir, "pod5")
            output_dir.mkdirs()
            it.pod5.each { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_pod5.map {
            if (it.get("dorado_job_bytes") || it.get("dorado_jobs")) {
                [*:it, pod5_unchunked: it["pod5"], pod5: chunk_files(it["pod5"], it.get("dorado_jobs"), it.get("dorado_job_bytes"))]
            } else {
                it
            }
        }
    .set { ch_pod5_chunked }
    call_process(DORADO_DOWNLOAD,
        ch_pod5_chunked,
        ["dorado_model", "dorado_download_args"],
        [id: { it.dorado_model }],
        ["dorado_model_dir"]) { meta -> [meta.dorado_model] }
    .set { ch_dorado_model }
    ch_dorado_model.branch {
            yes: !(it.containsKey("duplex") && !it["duplex"])
            no: true
        }
    .set { ch_do_duplex }
    call_process(DORADO_DOWNLOAD2,
        ch_do_duplex.yes,
        ["dorado_duplex_model", "dorado_download_args"],
        [id: { it.dorado_duplex_model }],
        ["dorado_duplex_model_dir"]) { meta -> [meta.dorado_duplex_model] }
    .set { ch_dorado_duplex_model }
    map_call_process(DORADO_DUPLEX,
        ch_dorado_duplex_model,
        ["pod5", "dorado_model_dir", "dorado_duplex_model_dir", "dorado_duplex_args"],
        [id: { pod5, meta -> exemplar(pod5).baseName }],
        "pod5",
        ["bam"],
        "_DORADO_DUPLEX") { pod5, meta -> [pod5, meta.dorado_model_dir, meta.dorado_duplex_model_dir] }
    .set { ch_dorado_duplex }
    map_call_process(DORADO_BASECALLER,
        ch_do_duplex.no,
        ["pod5", "dorado_model_dir", "dorado_basecaller_args"],
        [id: { pod5, meta -> exemplar(pod5).baseName }],
        "pod5",
        ["bam"],
        "_DORADO_BASECALLER") { pod5, meta -> [pod5, meta.dorado_model_dir] }
    .set { ch_dorado_basecaller }
    ch_dorado_basecaller.mix(ch_dorado_duplex)
    .set { ch_basecalled }
    ch_basecalled.subscribe {
        if (it.get("publish_bam")) {
            def output_dir = file_in_dir(it.output_run_dir, "bam")
            output_dir.mkdirs()
            it.bam.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_basecalled.mix(ch_input_type.bam.map { [*:it, bam: it.bam_input] } )
    .set { ch_bam }
    // use samtools to convert sam to fastq.gz
    map_call_process(SAMTOOLS_FASTQ,
        ch_bam,
        ["bam", "samtools_fastq_args"],
        [id: { bam, meta -> bam.baseName }],
        "bam",
        ["fastq"],
        "_SAMTOOLS_FASTQ") { bam, meta -> [bam] }
    .set { ch_converted_to_fastq }
    // publish fastq
    ch_converted_to_fastq.subscribe {
        if (it.get("publish_fastq")) {
            def output_dir = file_in_dir(it.output_run_dir, "fastq")
            output_dir.mkdirs()
            it.fastq.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_converted_to_fastq.mix(ch_input_type.fastq.map { [*:it, fastq: it.fastq_input] })
    .set { ch_fastq }
    ch_fastq.branch {
        yes: it.getOrDefault("align", true)
        no: true
    }
    .set { ch_to_align }
    // // GraphAligner -t 1 -x dbg -g barcode.gfa -f in.fastq.gz -a out.gaf
    with_keys(ch_to_align.yes, [graphaligner_args: { it.graphaligner_grouping_args }]) {
        map_call_process(GRAPHALIGNER_GROUPING,
            it,
            ["fastq", "gfa_grouping", "graphaligner_args"],
            [id: { fastq, meta -> fastq.name.replaceFirst(/\.fastq\.gz$/, "") }],
            "fastq",
            ["gaf_grouping"],
            "_GRAPHALIGNER_GROUPING") { fastq, meta -> [fastq, meta.gfa_grouping] }
    }
    .map {
        // bam, gaf -> [(bam, gaf), ...]
        [*:it, bam_and_gaf: [it.bam, it.gaf_grouping].transpose()]
    }
    .set { ch_graphaligner_grouping }
    // bin/join_gaf.py a.bam a.gaf out.arrow
    with_keys(ch_graphaligner_grouping, [join_gaf_args: { it.join_gaf_grouping_args }]) {
        map_call_process(JOIN_GAF_GROUPING,
            it,
            ["bam_and_gaf", "join_gaf_args", "tabular_format"],
            [
                id: { bam_and_gaf, meta -> bam_and_gaf[0].baseName },
                input_format: { bam_and_gaf, meta -> "bam" },
                output_format: { bam_and_gaf, meta -> meta.tabular_format }
            ],
            "bam_and_gaf",
            ["join_gaf_grouping_output"],
            "_JOIN_GAF_GROUPING") { bam_and_gaf, meta -> [bam_and_gaf[0], bam_and_gaf[1]] }
    }
    .set { ch_join_gaf_grouping }
    // bin/prepare_reads.py --gfa barcode.gfa in.arrow out.arrow -x UNS9,BC:UPSTREAM,BC:JUNCTION,BC:T7_TERM,BC:SPACER2
    map_call_process(PREPARE_READS,
        ch_join_gaf_grouping,
        ["join_gaf_grouping_output", "gfa_grouping", "prepare_reads_args", "tabular_format"],
        [
            id: { join_gaf_grouping_output, meta -> join_gaf_grouping_output.baseName },
            input_format: { join_gaf_grouping_output, meta -> meta.tabular_format },
            output_format: { join_gaf_grouping_output, meta -> meta.tabular_format }
        ],
        "join_gaf_grouping_output",
        ["prepare_reads_output"],
        "_PREPARE_READS") { join_gaf_grouping_output, meta -> [join_gaf_grouping_output, meta.gfa_grouping] }
    .set { ch_did_prepare_reads }
    ch_did_prepare_reads.subscribe {
        if (it.getOrDefault("publish_prepare_reads", true)) {
            def output_dir = file_in_dir(it.output_run_dir, "prepare_reads_output")
            output_dir.mkdirs()
            it.prepare_reads_output.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_did_prepare_reads.mix(ch_input_type.prepare_reads.map { [*:it, prepare_reads_output: it.prepare_reads_input] })
    .set { ch_prepare_reads }
    // set consensus group size based on file size?
    ch_prepare_reads.map {
        [*:it, consensus_groups: (0..<it.consensus_jobs).toList()]
    }
    .set { ch_consensus_groups }
    // bin/consensus.py "*.arrow" --output out8.arrow --fasta out8.fasta --group 8/100 --min-depth 5 --no-phred-output (??) --method-spoa
    map_call_process(CONSENSUS,
        ch_consensus_groups,
        ["consensus_groups", "prepare_reads_output", "consensus_args", "tabular_format"],
        [
            id: { consensus_group, meta -> "consensus-${consensus_group}-of-${meta.consensus_groups.size}" },
            group: { consensus_group, meta -> "${consensus_group}/${meta.consensus_groups.size}" },
            input_format: { consensus_group, meta -> meta.tabular_format },
            output_format: { consensus_group, meta -> meta.tabular_format }
        ],
        "consensus_groups",
        ["consensus_tabular", "consensus_fasta"],
        "_CONSENSUS") { consensus_group, meta -> [meta.prepare_reads_output] }
    .set { ch_did_consensus }
    ch_did_consensus.subscribe {
        if (it.get("publish_consensus")) {
            def output_dir = file_in_dir(it.output_run_dir, "consensus")
            output_dir.mkdirs()
            it.consensus_tabular.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
            it.consensus_fasta.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_input_type.consensus.map {
        def fasta_map = it.consensus_fasta_input.collectEntries { f -> [(f.baseName): f] }
        def reordered_fasta = it.consensus_tabular_input.collect { f ->
            if (!fasta_map.containsKey(f.baseName)) {
                throw new Exception("could not find FASTA corresponding to tabular consensus input: ${f.name}")
            }
            fasta_map[f.baseName]
        }
        [*:it, consensus_tabular: it.consensus_tabular_input, consensus_fasta: reordered_fasta]
    }
    .set { ch_input_consensus }
    ch_did_consensus.mix(ch_input_consensus)
    .set { ch_consensus }
    // regroup (minion: 100 consensus groups -> 10 new groups)
    // need to synchronize fasta regrouping with arrow regrouping
    ch_consensus.map {
        def align_chunks = [it.consensus_tabular, it.consensus_fasta].transpose().collate(it.consensus_jobs_per_align_job.toInteger())
        def consensus_tabular = align_chunks.collect { chunks -> chunks.collect { files -> files[0] } }
        def consensus_fasta = align_chunks.collect { chunks -> chunks.collect { files -> files[1] } }
        [*:it, consensus_tabular: consensus_tabular, consensus_fasta: consensus_fasta]
    }
    .set { ch_consensus_rechunked }
    // GraphAligner -t 1 -x dbg out3.fasta out3.gaf
    with_keys(ch_consensus_rechunked, [graphaligner_args: { it.graphaligner_variants_args }]) {
        map_call_process(GRAPHALIGNER_VARIANTS,
            it,
            ["consensus_fasta", "gfa_variants", "graphaligner_args"],
            [id: { fasta, meta -> exemplar(fasta).baseName }],
            "consensus_fasta",
            ["gaf_variants"],
            "_GRAPHALIGNER_VARIANTS") { fasta, meta -> [fasta, meta.gfa_variants] }
    }
    .map {
        // bam, gaf -> [(bam, gaf), ...]
        [*:it, consensus_and_gaf: [it.consensus_tabular, it.gaf_variants].transpose()]
    }
    .set { ch_graphaligner_variants }
    // bin/join_gaf.py --gaf out3.gaf --reads-prefix read_group_ --gaf-prefix consensus_
    with_keys(ch_graphaligner_variants, [join_gaf_args: { it.join_gaf_variants_args }]) {
        map_call_process(JOIN_GAF_VARIANTS,
            it,
            ["consensus_and_gaf", "join_gaf_args", "tabular_format"],
            [
                id: { consensus_and_gaf, meta -> exemplar(consensus_and_gaf)[0].baseName },
                input_format: { consensus_and_gaf, meta -> meta.tabular_format },
                output_format: { consensus_and_gaf, meta -> meta.tabular_format }
            ],
            "consensus_and_gaf",
            ["join_gaf_variants_output"],
            "_JOIN_GAF_VARIANTS") { consensus_and_gaf, meta -> [consensus_and_gaf[0], consensus_and_gaf[1]] }
    }
    .set { ch_join_gaf_variants }
    // bin/realign.py --gfa full.gfa out3.arrow out3_realigned.arrow
    map_call_process(REALIGN,
        ch_join_gaf_variants,
        ["join_gaf_variants_output", "gfa_variants", "realign_args", "tabular_format"],
        [
            id: { join_gaf_variants_output, meta -> join_gaf_variants_output.baseName },
            input_format: { join_gaf_variants_output, meta -> meta.tabular_format },
            output_format: { join_gaf_variants_output, meta -> meta.tabular_format }
        ],
        "join_gaf_variants_output",
        ["realign_output"],
        "_REALIGN") { join_gaf_variants_output, meta -> [join_gaf_variants_output, meta.gfa_variants] }
    .set { ch_did_realign }
    ch_did_realign.subscribe {
        if (it.get("publish_realign")) {
            def output_dir = file_in_dir(it.output_run_dir, "realign")
            output_dir.mkdirs()
            it.realign_output.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_did_realign.mix(ch_input_type.realign.map { [*:it, realign_output: it.realign_input] })
    .set { ch_realign }
    // bin/join_gaf.py "*.arrow" combined.arrow
    // bin/extract_segments.py --path-col consensus_path --cigar-col realign_cg --gfa full.gfa realigned3.arrow realigned3_extracted.arrow
    map_call_process(EXTRACT_SEGMENTS,
        ch_realign,
        ["realign_output", "gfa_variants", "extract_segments_args", "tabular_format"],
        [
            id: { realign_output, meta -> realign_output.baseName },
            input_format: { realign_output, meta -> meta.tabular_format },
            output_format: { realign_output, meta -> meta.tabular_format }
        ],
        "realign_output",
        ["extract_segments_output"],
        "_EXTRACT_SEGMENTS") { realign_output, meta -> [realign_output, meta.gfa_variants] }
        .set { ch_extract_segments }
    // publish tabular output
    ch_extract_segments.subscribe {
        if (it.getOrDefault("publish_extract_segments", true)) {
            def output_dir = file_in_dir(it.output_run_dir, "extract_segments")
            output_dir.mkdirs()
            it.extract_segments_output.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_extract_segments.mix(ch_to_align.no)
        .set { samples }

    samples.view()

    emit:
    samples
}

workflow MAIN {
    samples_in = get_samples(params, [:], true)
    samples_in = find_inputs(samples_in, params.root, ["gfa_grouping", "gfa_variants", "gfa"])
    samples_in = glob_inputs(samples_in, params.root, GLOBBED_INPUTS)
    samples = NANOPORE_FISH(samples_in)

    emit:
    samples
}

workflow {
    MAIN()
}
