import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import static functions.*

include { POD5_MERGE; POD5_MERGE as POD5_MERGE2; POD5_VIEW; POD5_FILTER; SPLIT_READ_IDS } from '../../modules/pod5.nf'

workflow NANOPORE_FISH {
    take:
    samples_in

    main:
    samples_in
        .branch {
            // chunk pod5 files UNLESS:
            // map contains the key pod5_chunk and it is false
            // OR both pod5_chunk_bytes and pod5_chunk_files are falsy
            yes: ! ((it.containsKey("pod5_chunk") && !it["pod5_chunk"])
                                        || (!it.get("pod5_chunk_bytes")
                                            && !it.get("pod5_chunk_files")))
            no: true
        }
        . set { ch_do_pod5_chunk }
    ch_do_pod5_chunk.yes.map {
            if (it.get("pod5_input")) {
                it["pod5_input_chunked"] = chunk_files(it["pod5_input"], it.get("pod5_chunk_bytes"), it.get("pod5_chunk_files"))
            }
            it
        }
        . set { ch_pod5_chunked }
    ch_pod5_chunked.set { samples }
    map_call_process(POD5_MERGE,
                     ch_pod5_chunked,
                     ["pod5_input_chunked", "pod5_merge_args"],
                     [id: { pod5, meta -> pod5[0].baseName }],
                     "pod5_input_chunked",
                     ["pod5_to_split"],
                     "_POD5_MERGE") { pod5, meta -> [pod5] }
        . set { ch_pod5_merged }
    ch_pod5_merged.mix(ch_do_pod5_chunk.no.map {
                [*:it, pod5_to_split: it.get("pod5_input")]
            } )
        .branch {
            yes: !(it.containsKey("pod5_split") && !it["pod5_split"])
            no: true
        }
        .set { ch_do_pod5_split }
    ch_do_pod5_split.yes.map {
            def pod5_split_by = it.get("pod5_split_by") ?: ["channel"]
            [*:it,
             pod5_split_by: pod5_split_by,
             split_read_ids_args: "-F \"${pod5_split_by.join(',')}\" ${it.get('split_read_ids_args') ?: ''}",
             pod5_view_args: "--include \"${['read_id', *pod5_split_by].join(',')}\" ${it.get('pod5_view_args') ?: ''}"]
        }
        .set { ch_to_pod5_split }
    call_process(POD5_VIEW,
        ch_to_pod5_split,
        ["pod5_to_split", "pod5_view_args"],
        [id: { it.pod5_to_split[0].baseName }],
        ["pod5_read_list"]) { meta -> [meta.pod5_to_split] }
        .set { ch_pod5_read_lists }
    call_process(SPLIT_READ_IDS,
        ch_pod5_read_lists,
        ["pod5_read_list", "split_read_ids_args"],
        [id: { it.pod5_read_list.baseName }],
        ["pod5_split_read_lists"]) { meta -> [meta.pod5_read_list] }
        .set { ch_pod5_split_read_lists }
    map_call_process(POD5_FILTER,
        ch_pod5_split_read_lists.map { edit_key(it, "pod5_split_read_lists") { l -> l[0..10] } },
        // ch_pod5_split_read_lists,
        ["pod5_to_split", "pod5_filter_args"],
        [id: { read_list, meta -> read_list.baseName }],
        "pod5_split_read_lists",
        ["pod5_split"],
        "_POD5_FILTER") { read_list, meta -> [meta.pod5_to_split, read_list] }
        .set { samples }

    // dorado simplex+duplex model download
    // dorado duplex basecalling
    // publish bams
    // filter dx:0/1 from bam -> fastq.gz
    // porechop_abi
    // minigraph barcode_v1.gfa
    // parse gaf, output mapping of barcode to list of read ids (using pyarrow dump)
    // index fastq
    // download medaka model
    // consensus.py
    // identify_variants.py
    // merge parquets

        // | map {
        //     if (it.get("filtered_reads")) {
        //         rename_key(rename_key(it, "reads", "unfiltered_reads"), "filtered_reads", "reads")
        //     } else {
        //         it
        //     }
        // }
        // | call_MINIMAP2_ALIGN
        // | call_SAMTOOLS_INDEX
        // | map {
        //     it.output_run_dir.mkdirs()
        //     ["reference", "reads", "bam", "bam_index"].each { k ->
        //         (it[k] instanceof Collection ? it[k] : [it[k]]).each { src ->
        //             Files.copy(src, file_in_dir(it.output_run_dir, src.name),
        //                        StandardCopyOption.REPLACE_EXISTING)
        //         }
        //     }
        //     it
        // }
        // | set { samples }

    samples.view()

    emit:
    samples

}

workflow MAIN {
    // download_data(params) // TODO: remove?
    samples_preglob = get_samples(params, [:], true)
    glob_inputs(samples_preglob, params.root, ["fastq_input", "pod5_input"])
        | NANOPORE_FISH
        | set { samples }

    emit:
    samples
}

workflow {
    MAIN()
}
