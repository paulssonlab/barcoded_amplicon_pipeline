import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import static functions.*

include { PREPARE_REFERENCES } from '../prepare.nf'

include { call_MINIMAP2_INDEX;
          call_MINIMAP2_ALIGN } from '../../modules/minimap2.nf'

include { call_CHOPPER } from '../../modules/chopper.nf'

include { call_SAMTOOLS_INDEX } from '../../modules/samtools.nf'

include { call_CALL_VARIANTS;
          call_FILTER_VARIANTS;
          call_INDEX_VARIANTS;
          call_GET_CONSENSUS } from '../../modules/variants.nf'

include { call_EXTRACT_CONSENSUS } from '../../modules/fastas.nf'

workflow NANOPORE {
    take:
    samples_in

    main:
    samples_in
        | call_MINIMAP2_INDEX
        | call_CHOPPER
        | map {
            if (it.get("filtered_reads")) {
                rename_key(rename_key(it, "reads", "unfiltered_reads"), "filtered_reads", "reads")
            } else {
                it
            }
        }
        | call_MINIMAP2_ALIGN
        | call_SAMTOOLS_INDEX
        | map {
            it.output_run_dir.mkdirs()
            ["reference", "reads", "bam", "bam_index"].each { k ->
                (it[k] instanceof Collection ? it[k] : [it[k]]).each { src ->
                    Files.copy(src, file_in_dir(it.output_run_dir, src.name),
                               StandardCopyOption.REPLACE_EXISTING)
                }
            }
            it
        }
        // | call_CALL_VARIANTS
        // | call_FILTER_VARIANTS
        // | call_INDEX_VARIANTS
        // | call_GET_CONSENSUS
        // | call_EXTRACT_CONSENSUS
        | set { samples }

    samples.view()

    emit:
    samples

}

workflow MAIN {
    download_data(params)
    samples_preglob = get_samples(params, [:], true) {
        if (it.get("ignore_references")) {
            it.references = ""
        }
    }
    glob_inputs(samples_preglob, params.data_dir, ["fastq", "fast5", "pod5"])
        | PREPARE_REFERENCES
        | map { rename_key(it, "fastq", "reads") }
        | NANOPORE
        | set { samples }

    emit:
    samples
}

workflow {
    MAIN()
}
