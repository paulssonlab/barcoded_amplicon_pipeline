import static functions.*

include { PREPARE_SAMPLE_SHEET;
          PREPARE_READS;
          PREPARE_REFERENCES } from '../prepare.nf'

include { call_BOWTIE2_BUILD;
          call_BOWTIE2_ALIGN } from '../../modules/bowtie2.nf'

include { call_SAMTOOLS_INDEX } from '../../modules/samtools.nf'
include { call_CALL_VARIANTS;
          call_FILTER_VARIANTS;
          call_INDEX_VARIANTS;
          call_GET_CONSENSUS } from '../../modules/variants.nf'

include { call_EXTRACT_CONSENSUS } from '../../modules/fastas.nf'

workflow ILLUMINA_WHOLE_PLASMID {
    take:
    samples_in

    main:
    PREPARE_REFERENCES(samples_in)
        | call_BOWTIE2_BUILD
        | call_BOWTIE2_ALIGN
        | call_SAMTOOLS_INDEX
        | call_CALL_VARIANTS
        | call_FILTER_VARIANTS
        | call_INDEX_VARIANTS
        | call_GET_CONSENSUS
        | call_EXTRACT_CONSENSUS
        | set { samples }

    emit:
    samples

}

workflow MAIN {
    PREPARE_SAMPLE_SHEET()
        // map reads_prefix to reads_path
        .map { [reads_path: "${it.reads_prefix}.fastq", *:it] }
        | PREPARE_READS
        | ILLUMINA_WHOLE_PLASMID
        | set { samples }

    emit:
    samples
}

workflow {
    MAIN()
}
