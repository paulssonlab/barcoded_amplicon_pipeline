import static functions.*

include { PREPARE_SAMPLE_SHEET;
          PREPARE_READS;
          PREPARE_REFERENCES } from '../prepare.nf'

include { call_BOWTIE2_BUILD;
          call_BOWTIE2_ALIGN } from '../../modules/bowtie2.nf'

include { call_SAMTOOLS_SORT;
          call_SAMTOOLS_INDEX } from '../../modules/samtools.nf'

workflow ILLUMINA_WHOLE_PLASMID {
    take:
    samples_in

    main:
    PREPARE_REFERENCES(samples_in)
        | call_BOWTIE2_BUILD
        | call_BOWTIE2_ALIGN
        | view
    // call_BOWTIE2_BUILD(PREPARE_REFERENCES.out)
    //     .view()

    // BOWTIE2_BUILD(PREPARE_REFERENCES.out.references)
    // MERGE_INDEXES(PREPARE_REFERENCES.out.samples, BOWTIE2_BUILD.out.index)
    // MERGE_INDEXES.out.samples.view()

    // call_BOWTIE2(MERGE_INDEXES.out.samples)
    //     .set { ch_mapped }

    // call_SAMTOOLS_SORT(ch_mapped)
    //     .set { ch_sorted }

    // call_SAMTOOLS_INDEX(ch_sorted)
    //     .set { ch_indexed }

    // ch_indexed
    //     .view()
        // .set { ch_sorted }

    // CALL_VARIANTS
    // FILTER_VARIANTS
    // INDEX_VARIANTS (?)
    // GET_CONSENSUS
    // EXTRACT_CONSENSUS

    // TODO: fix /tmp/.../tmp issue

    //emit:
    //samples

}

workflow MAIN {
    PREPARE_SAMPLE_SHEET().samples
        // map reads_prefix to reads_path
        .map { [reads_path: "${it.reads_prefix}.fastq", *:it] }
        | PREPARE_READS
        | ILLUMINA_WHOLE_PLASMID
}

workflow {
    MAIN()
}
