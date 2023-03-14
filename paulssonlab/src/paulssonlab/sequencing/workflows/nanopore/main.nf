import static functions.*

include { PREPARE_REFERENCES } from '../prepare.nf'

include { call_MINIMAP2_INDEX;
          call_MINIMAP2_ALIGN } from '../../modules/minimap2.nf'

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
        // | call_MINIMAP2_INDEX
        // | call_MINIMAP2_ALIGN
        // | call_SAMTOOLS_INDEX
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
        | map { renameKey(it, "fastq", "reads") }
        | NANOPORE
        | set { samples }

    emit:
    samples
}

workflow {
    MAIN()
}
