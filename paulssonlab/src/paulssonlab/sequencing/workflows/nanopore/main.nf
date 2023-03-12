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
        | call_MINIMAP2_INDEX
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
    // println glob("*.fastq.gz", "/tmp/paulssonlab-sequencing/230308_repressilator_debugging/data")
    // download_data(params)
    glob_inputs(get_samples(params), params.data_dir, ["fastq", "fast5", "pod5"])
        // | PREPARE_REFERENCES
        | view
        // | NANOPORE
        // | set { samples }

    // emit:
    // samples
}

workflow {
    MAIN()
}
