import static functions.*

include { POD5_MERGE; POD5_MERGE as POD5_MERGE2; POD5_VIEW_AND_SUBSET } from '../../modules/pod5.nf'
include { DORADO_DOWNLOAD;
          DORADO_DOWNLOAD as DORADO_DOWNLOAD2;
          DORADO_BASECALLER;
          DORADO_DUPLEX;
          DORADO_DUPLEX_WITH_PAIRS } from '../../modules/dorado.nf'
include { SAMTOOLS_FASTQ;
          SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_DUPLEX;
          SAMTOOLS_IMPORT;
          SAMTOOLS_MERGE } from '../../modules/samtools.nf'
include { GRAPHALIGNER as GRAPHALIGNER_GROUPING;
          GRAPHALIGNER as GRAPHALIGNER_GROUPING_DUPLEX;
          GRAPHALIGNER as GRAPHALIGNER_VARIANTS } from '../../modules/graphaligner.nf'
include { FIND_DUPLEX_PAIRS;
          JOIN_GAF as JOIN_GAF_GROUPING;
          JOIN_GAF as JOIN_GAF_VARIANTS;
          PREPARE_READS;
          CONSENSUS;
          PREPARE_CONSENSUS;
          CONSENSUS_PREPARED;
          REALIGN;
          EXTRACT_SEGMENTS } from '../../modules/scripts.nf'
include { CAT as CAT_DUPLEX_GAF;
          CAT as CAT_DUPLEX_FASTQ } from '../../modules/util.nf'

def GLOBBED_INPUTS = ["pod5_input", "bam_input", "fastq_input", "prepare_reads_input", "prepare_consensus_input", "consensus_tabular_input", "consensus_fasta_input", "realign_input"]
def REQUIRED_INPUTS = ["bam_input", "fastq_input", "prepare_reads_input", "prepare_consensus_input", "consensus_tabular_input", "realign_input"]

workflow NANOPORE_FISH {
    take:
    samples_in

    main:
    def DEFAULT_ARGS = [
        basecall: true,
        duplex: true,
        use_dorado_duplex_pairing: false,
        pod5_chunk: true,
        publish_pod5: true,
        publish_bam: true,
        publish_fastq: false,
        publish_prepare_consensus: false,
        publish_realign: false,
        publish_extract_segments: true,
        output: "extract_segments", // possible values: "pod5", "basecaller", "consensus", "extract_segments"
        samtools_merge_args: "-c",
        tabular_format: "arrow",
        graphaligner_args: "-x dbg",
        find_duplex_pairs_args: "-x UNS9,BC:UPSTREAM,BC:JUNCTION,BC:T7_TERM,BC:SPACER2",
        prepare_reads_args: "-x UNS9,BC:UPSTREAM,BC:JUNCTION,BC:T7_TERM,BC:SPACER2",
        prepare_consensus: true,
        consensus_args: "--method spoa --no-phred-output --min-depth 3",
        consensus_jobs: 200,
        consensus_jobs_per_align_job: 1,
        join_gaf_variants_args: "--rename-col path grouping_path --rename-gaf-col path variants_path",
        //extract_segments_args: "--path-col variants_path --cigar-col realign_cg",
    ]
    samples_in
    .map {
        if (it.get("basecall") && !it.get("pod5_input")) {
            throw new Exception("basecalling requires pod5 input")
        }
        if (it.get("bam_input") && it.get("fastq_input")) {
            throw new Exception("cannot specify both bam_input and fastq_input")
        }
        if ((!it.get("consensus_tabular_input") && it.get("consensus_fasta_input")) || (it.get("consensus_tabular_input") && !it.get("consensus_fasta_input"))) {
            throw new Exception("both consensus_tabular_input and consensus_fasta_input must be given to start pipeline from consensus step")
        }
        if (!["basecall", *REQUIRED_INPUTS].collect { k -> it.get(k) }.any()) {
            throw new Exception("missing any pipeline input")
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
        if (it["align"]) {
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
        prepare_consensus: it.get("prepare_consensus_input")
        prepare_reads: it.get("prepare_reads_input")
        bam: it.get("bam_input")
        fastq: it.get("fastq_input")
    }
    .set { ch_input_type }
    ch_input_type.pod5.branch {
        // chunk pod5 files UNLESS:
        // map contains the key pod5_chunk and it is false
        // OR both pod5_chunk_bytes and pod5_chunk_files are falsy
        yes: !(!it["pod5_chunk"]
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
        if (it["publish_pod5"]) {
            def output_dir = file_in_dir(it.output_run_dir, "pod5")
            output_dir.mkdirs()
            it.pod5.each { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_pod5.branch {
        no: it["output"] == "pod5"
        yes: true
    }
    .set { ch_process_pod5 }
    ch_process_pod5.yes.map {
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
            yes: it["duplex"]
            no: true
        }
    .set { ch_do_download_duplex }
    call_process(DORADO_DOWNLOAD2,
        ch_do_download_duplex.yes,
        ["dorado_duplex_model", "dorado_download_args"],
        [id: { it.dorado_duplex_model }],
        ["dorado_duplex_model_dir"]) { meta -> [meta.dorado_duplex_model] }
    .set { ch_did_download_duplex }
    ch_did_download_duplex.branch {
        yes: it["duplex"] && it["use_dorado_duplex_pairing"]
        no: true
    }
    .set { ch_do_dorado_duplex }
    map_call_process(DORADO_DUPLEX,
        ch_do_dorado_duplex.yes,
        ["pod5", "dorado_model_dir", "dorado_duplex_model_dir", "dorado_duplex_args"],
        [id: { pod5, meta -> exemplar(pod5).baseName }],
        "pod5",
        ["bam"],
        "_DORADO_DUPLEX") { pod5, meta -> [pod5, meta.dorado_model_dir, meta.dorado_duplex_model_dir] }
    .set { ch_did_dorado_duplex }
    ch_do_download_duplex.no.mix(ch_do_dorado_duplex.no)
    .set { ch_do_dorado_simplex }
    map_call_process(DORADO_BASECALLER,
        ch_do_dorado_simplex,
        ["pod5", "dorado_model_dir", "dorado_basecaller_args"],
        [id: { pod5, meta -> exemplar(pod5).baseName }],
        "pod5",
        ["bam"],
        "_DORADO_BASECALLER") { pod5, meta -> [pod5, meta.dorado_model_dir] }
    .set { ch_did_dorado_simplex }
    ch_did_dorado_simplex.mix(ch_did_dorado_duplex)
    .set { ch_basecalled }
    // publish bam
    ch_basecalled.subscribe {
        if (it["publish_bam"] && it["use_dorado_duplex_pairing"]) {
            def output_dir = file_in_dir(it.output_run_dir, "bam")
            output_dir.mkdirs()
            it.bam.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_basecalled.mix(ch_input_type.bam.map { [*:it, bam: it.bam_input, basecall: false] } )
    .set { ch_bam }
    // use samtools to convert sam to fastq.gz
    map_call_process(SAMTOOLS_FASTQ,
        ch_bam,
        ["bam", "samtools_fastq_args"],
        [id: { bam, meta -> bam.baseName }],
        "bam",
        ["fastq"],
        "_SAMTOOLS_FASTQ") { bam, meta -> [bam] }
    .set { ch_converted_bam_to_fastq }
    // publish fastq
    ch_converted_bam_to_fastq.subscribe {
        if (it["publish_fastq"] && it["use_dorado_duplex_pairing"]) {
            def output_dir = file_in_dir(it.output_run_dir, "fastq")
            output_dir.mkdirs()
            it.fastq.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_input_type.fastq.map { [*:it, fastq: it.fastq_input, basecall: false] }
    .set { ch_fastq_input }
    // convert input fastq to bam
    map_call_process(SAMTOOLS_IMPORT,
        ch_fastq_input,
        ["fastq", "samtools_import_args"],
        [id: { fastq, meta -> fastq.baseName }],
        "fastq",
        ["bam"],
        "_SAMTOOLS_IMPORT") { fastq, meta -> [fastq] }
    .set { ch_converted_fastq_to_bam }
    ch_converted_bam_to_fastq.mix(ch_converted_fastq_to_bam)
    .set { ch_bam_and_fastq }
    ch_bam_and_fastq.branch {
        no: it["output"] == "basecaller"
        yes: true
    }
    .set { ch_process_basecalled_reads }
    // GraphAligner -t 1 -x dbg -g barcode.gfa -f in.fastq.gz -a out.gaf
    with_keys(ch_process_basecalled_reads.yes, [graphaligner_args: { it.graphaligner_grouping_args }]) {
        map_call_process(GRAPHALIGNER_GROUPING,
            it,
            ["fastq", "gfa_grouping", "graphaligner_args"],
            [id: { fastq, meta -> fastq.name.replaceFirst(/\.fastq\.gz$/, "") }],
            "fastq",
            ["gaf_grouping"],
            "_GRAPHALIGNER_GROUPING") { fastq, meta -> [fastq, meta.gfa_grouping] }
    }
    .set { ch_did_graphaligner_grouping }
    ch_did_graphaligner_grouping.map {
        // bam, gaf -> [(bam, gaf), ...]
        [*:it, bam_and_gaf: [it.bam, it.gaf_grouping].transpose()]
    }
    .set { ch_bam_and_gaf }
    ch_bam_and_gaf.branch {
        yes: it["basecall"] && it["duplex"] && !it["use_dorado_duplex_pairing"]
        no: true
    }
    .set { ch_do_nondorado_duplex_pairing }
    // start of use_dorado_duplex_pairing = false
    // find_duplex_pairs.py --gfa ../references/barcode.gfa --gaf channel-1_merged.gaf -x UNS9,BC:UPSTREAM,BC:JUNCTION,BC:T7_TERM,BC:SPACER2 channel-1_merged.bam channel-1_merged_pairs.csv
    map_call_process(FIND_DUPLEX_PAIRS,
        ch_do_nondorado_duplex_pairing.yes,
        ["bam_and_gaf", "gfa_grouping", "find_duplex_pairs_args"],
        [id: { bam_and_gaf, meta -> bam_and_gaf[0].baseName }],
        "bam_and_gaf",
        ["duplex_pairs"],
        "_FIND_DUPLEX_PAIRS") { bam_and_gaf, meta -> [meta.gfa_grouping, bam_and_gaf[0], bam_and_gaf[1]] }
    .set { ch_did_find_duplex_pairs }
    ch_did_find_duplex_pairs.map {
        [*:it, pod5_and_pairs: [it.pod5, it.duplex_pairs].transpose()]
    }
    .set { ch_pod5_and_pairs }
    // dorado duplex --pairs channel-1_merged_pairs.csv simplex_model pod5
    map_call_process(DORADO_DUPLEX_WITH_PAIRS,
        ch_pod5_and_pairs,
        ["pod5_and_pairs", "dorado_model_dir", "dorado_duplex_model_dir", "dorado_duplex_args"],
        [id: { pod5_and_pairs, meta -> pod5_and_pairs[0].baseName }],
        "pod5_and_pairs",
        ["bam_duplex"],
        "_DORADO_DUPLEX_WITH_PAIRS") { pod5_and_pairs, meta -> [pod5_and_pairs[0], pod5_and_pairs[1], meta.dorado_model_dir, meta.dorado_duplex_model_dir] }
    .set { ch_did_dorado_duplex }
    // samtools merge -@ ${task.cpus} -c bam/channel-1_merged.bam bam/channel-1_merged_duplex.bam channel-1_merged.bam
    ch_did_dorado_duplex.map {
        [
            *:it,
            bam_simplex: it.bam,
            bam_to_combine: [it.bam_duplex, it.bam].transpose()
        ]
    }
    .set { ch_bam_to_combine }
    map_call_process(SAMTOOLS_MERGE,
        ch_bam_to_combine,
        ["bam_to_combine", "samtools_merge_args"],
        [id: { bam_to_combine, meta -> bam_to_combine[1].baseName }],
        "bam_to_combine",
        ["bam"],
        "_SAMTOOLS_MERGE") { bam_to_combine, meta -> [bam_to_combine] }
    .set { ch_bam_combined }
    // publish bam
    ch_bam_combined.subscribe {
        if (it["publish_bam"]) {
            def output_dir = file_in_dir(it.output_run_dir, "bam")
            output_dir.mkdirs()
            it.bam.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    // use samtools to convert sam to fastq.gz
    map_call_process(SAMTOOLS_FASTQ_DUPLEX,
        ch_bam_combined,
        ["bam_duplex", "samtools_fastq_args"],
        [id: { bam, meta -> bam.baseName }],
        "bam_duplex",
        ["fastq_duplex"],
        "_SAMTOOLS_FASTQ_DUPLEX") { bam, meta -> [bam] }
    .set { ch_duplex_converted_to_fastq }
    ch_duplex_converted_to_fastq.map {
        [
            *:it,
            fastq_simplex: it.fastq,
            fastq_to_combine: [it.fastq_duplex, it.fastq].transpose()
        ]
    }
    .set { ch_fastq_to_combine }
    // cat fastq/channel-1_merged_duplex.fastq.gz fastq/channel-1_merged.fastq.gz > channel-1_merged.fastq.gz
    map_call_process(CAT_DUPLEX_FASTQ,
        ch_fastq_to_combine,
        ["fastq_to_combine"],
        [id: { fastq_to_combine, meta -> fastq_to_combine[1] }],
        "fastq_to_combine",
        ["fastq"],
        "_CAT_DUPLEX_FASTQ") { fastq_to_combine, meta -> [fastq_to_combine] }
    .set { ch_fastq_combined }
    // publish fastq
    ch_fastq_combined.subscribe {
        if (it["publish_fastq"]) {
            def output_dir = file_in_dir(it.output_run_dir, "fastq")
            output_dir.mkdirs()
            it.fastq.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    // GraphAligner
    with_keys(ch_fastq_combined, [graphaligner_args: { it.graphaligner_grouping_args }]) {
        map_call_process(GRAPHALIGNER_GROUPING_DUPLEX,
            it,
            ["fastq_duplex", "gfa_grouping", "graphaligner_args"],
            [id: { fastq, meta -> fastq.name.replaceFirst(/\.fastq\.gz$/, "") }],
            "fastq",
            ["gaf_grouping_duplex"],
            "_GRAPHALIGNER_GROUPING_DUPLEX") { fastq, meta -> [fastq, meta.gfa_grouping] }
    }
    .set { ch_did_graphaligner_grouping_duplex }
    ch_did_graphaligner_grouping_duplex.map {
        [
            *:it,
            gaf_grouping_simplex: it.gaf_grouping,
            gaf_to_combine: [it.gaf_grouping_duplex, it.gaf_grouping].transpose()
        ]
    }
    .set { ch_gaf_to_combine }
    // cat gaf/channel-1_merged_duplex.gaf gaf/channel-1_merged.gaf > channel-1_merged.gaf
    map_call_process(CAT_DUPLEX_GAF,
        ch_gaf_to_combine,
        ["gaf_to_combine"],
        [id: { gaf_to_combine, meta -> gaf_to_combine[1] }],
        "gaf_to_combine",
        ["gaf_grouping"],
        "_CAT_DUPLEX_GAF") { gaf_to_combine, meta -> [gaf_to_combine] }
    .set { ch_did_nondorado_duplex_pairing }
    // end of use_dorado_duplex_pairing = false
    ch_do_nondorado_duplex_pairing.no.mix(ch_did_nondorado_duplex_pairing)
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
    // publish prepared_reads
    ch_did_prepare_reads.subscribe {
        if (it["publish_prepare_reads"]) {
            def output_dir = file_in_dir(it.output_run_dir, "prepare_reads")
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
    ch_consensus_groups.branch {
        yes: it.get("prepare_consensus")
        no: true
    }
    .set { ch_should_prepare_consensus }
    map_call_process(PREPARE_CONSENSUS,
        ch_should_prepare_consensus.yes,
        ["consensus_groups", "prepare_reads_output", "prepare_consensus_args", "tabular_format"],
        [
            id: { consensus_group, meta -> "consensus-${consensus_group}-of-${meta.consensus_groups.size}" },
            group: { consensus_group, meta -> "${consensus_group}/${meta.consensus_groups.size}" },
            input_format: { consensus_group, meta -> meta.tabular_format },
            output_format: { consensus_group, meta -> meta.tabular_format }
        ],
        "consensus_groups",
        ["prepare_consensus_output"],
        "_PREPARE_CONSENSUS") { consensus_group, meta -> [meta.prepare_reads_output] }
    .set { ch_did_prepare_consensus }
    // publish prepare_consensus
    ch_did_prepare_consensus.subscribe {
        if (it["publish_prepare_consensus"]) {
            def output_dir = file_in_dir(it.output_run_dir, "prepare_consensus")
            output_dir.mkdirs()
            it.prepare_consensus_output.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_did_prepare_consensus.mix(ch_input_type.prepare_consensus.map { [*:it, prepare_consensus_output: it.prepare_consensus_input] })
    .set { ch_prepare_consensus }
    map_call_process(CONSENSUS_PREPARED,
        ch_prepare_consensus,
        ["prepare_consensus_output", "consensus_args", "tabular_format"],
        [
            id: { prepare_consensus_output, meta -> prepare_consensus_output.baseName },
            input_format: { prepare_consensus_output, meta -> meta.tabular_format },
            output_format: { prepare_consensus_output, meta -> meta.tabular_format }
        ],
        "prepare_consensus_output",
        ["consensus_tabular", "consensus_fasta"],
        "_CONSENSUS_PREPARED") { prepare_consensus_output, meta -> [prepare_consensus_output] }
    .set { ch_did_consensus_prepared }
    map_call_process(CONSENSUS,
        ch_should_prepare_consensus.no,
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
    .set { ch_did_consensus_unprepared }
    ch_did_consensus_unprepared.mix(ch_did_consensus_prepared)
    .set { ch_did_consensus }
    // publish consensus
    ch_did_consensus.subscribe {
        if (it["publish_consensus"]) {
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
    ch_did_consensus.branch {
        no: it["output"] == "consensus"
        yes: true
    }
    .set { ch_process_consensus }
    ch_process_consensus.yes.mix(ch_input_consensus)
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
    // publish realign
    ch_did_realign.subscribe {
        if (it["publish_realign"]) {
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
    // publish extract_segments
    ch_extract_segments.subscribe {
        if (it["publish_extract_segments"]) {
            def output_dir = file_in_dir(it.output_run_dir, "extract_segments")
            output_dir.mkdirs()
            it.extract_segments_output.collect { output_file ->
                output_file.toRealPath().mklink(file_in_dir(output_dir, output_file.name), overwrite: true)
            }
        }
    }
    ch_extract_segments.mix(ch_process_pod5.no,
                            ch_process_basecalled_reads.no,
                            ch_process_consensus.no)
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
