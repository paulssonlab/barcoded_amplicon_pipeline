import static functions.*

include { PREPARE_READS;
          PREPARE_REFERENCES } from '../prepare.nf'

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
    PREPARE_REFERENCES(samples_in)
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
    Channel.fromList(get_samples(params, [fastq: '${name}/*.fastq*', fast5: '${name}/*.fast5'], true))
    // Channel.fromList(get_samples(params, [fastq: '${name}/*.fastq*', fast5: '${name}/*.fast5'], true))
        // map reads_prefix to reads_path
        // .map { [reads_path: "${it.reads_prefix}.fastq", *:it] }
        | view
        // | PREPARE_READS
        // | NANOPORE
        // | set { samples }

    // emit:
    // samples
}

workflow {
    MAIN()
}
