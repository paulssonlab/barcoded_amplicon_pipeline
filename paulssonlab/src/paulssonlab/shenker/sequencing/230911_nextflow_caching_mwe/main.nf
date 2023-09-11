import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import static functions.*

include { POD5_MERGE; POD5_VIEW; POD5_FILTER; SPLIT_READ_IDS } from "${src}/sequencing/modules/pod5.nf"

workflow {
    main:
    def samples_in = Channel.of([pod5_input: [file("/home/jqs1/scratch/jqs1/sequencing/230818_bcd_rbses_subset/20230818_1343_1A_PAQ97606_f49ab41c/pod5_pass/PAQ97606_pass_f49ab41c_ad491743_999.pod5")]])
    map_call_process(POD5_MERGE,
                     samples_in,
                     ["pod5_input", "pod5_merge_args"],
                     [id: { pod5, meta -> pod5[0].baseName }],
                     "pod5_input",
                     ["pod5_to_split"],
                     "_POD5_MERGE") { pod5, meta -> [pod5] }
        .map {
            def pod5_split_by = it.get("pod5_split_by") ?: ["channel"]
            def pod5_split_chunks = it.get("pod5_split_chunks") ?: 400
            [*:it,
             pod5_view_args: "--include \"${['read_id', *pod5_split_by].join(',')}\" ${it.get('pod5_view_args') ?: ''}"]
        }
        .set { ch_to_pod5_split }
    call_process(POD5_VIEW,
        ch_to_pod5_split,
        ["pod5_to_split", "pod5_view_args"],
        [id: { it.pod5_to_split[0].baseName }],
        ["pod5_read_list"]) { meta -> [meta.pod5_to_split] }
        .set { samples }

    samples.view()

    emit:
    samples

}
