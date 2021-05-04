include { call_process } from '../../functions.nf'

include { PREPARE_SAMPLE_SHEET;
          PREPARE_READS;
          PREPARE_REFERENCES;
          MERGE_INDEXES } from '../prepare.nf'

include { BOWTIE2_BUILD;
          BOWTIE2 } from '../../modules/bowtie2.nf'

include { SAMTOOLS_SORT;
          SAMTOOLS_INDEX } from '../../modules/samtools.nf'

workflow ILLUMINA_WHOLE_PLASMID {
    take:
    samples_in

    main:
    PREPARE_REFERENCES(samples_in)
    // TODO: use call_process for BOWTIE2_BUILD so can build index with different arguments
    BOWTIE2_BUILD(PREPARE_REFERENCES.out.references)
    MERGE_INDEXES(PREPARE_REFERENCES.out.samples, BOWTIE2_BUILD.out.index)

    call_process(BOWTIE2,
                 MERGE_INDEXES.out.samples,
                 ["reads", "index", "bowtie2_args"],
                 [id: { it.reads.baseName }],
                 ["bam", "bowtie2_log"]) { [it.reads, it.index] }
        .set { ch_mapped }

    call_process(SAMTOOLS_SORT,
                 ch_mapped,
                 ["bam"],
                 [id: { it.bam.name }],
                 ["sorted_bam"]) { [it.bam] }
        .set { ch_sorted }

    // ch_sorted.map { it.sorted_bam.name }.view()
    call_process(SAMTOOLS_INDEX,
                 ch_sorted,
                 ["sorted_bam"],
                 [id: { it.sorted_bam.name }],
                 ["sorted_bam_index"]) { [it.sorted_bam] }
        .set { ch_indexed }

    ch_indexed
        .view()
        // .set { ch_sorted }

    // CALL_VARIANTS
    // FILTER_VARIANTS
    // INDEX_VARIANTS (?)
    // GET_CONSENSUS
    // EXTRACT_CONSENSUS

    // TODO: fix /tmp/.../tmp issue

    // emit:
    // samples

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

// include { GET_REGISTRY_SEQS;
//          SCP as DOWNLOAD_SAMPLES;
//          SCP as DOWNLOAD_READS;
//          ANY2FASTA;
//          MERGE_FASTAS;
//          BOWTIE2_BUILD;
//          BOWTIE2_INTERLEAVED;
//          SAMTOOLS_SORT;
//          SAMTOOLS_INDEX;
//          CALL_VARIANTS;
//          FILTER_VARIANTS;
//          INDEX_VARIANTS;
//          GET_CONSENSUS;
//          EXTRACT_CONSENSUS } from './modules.nf'

// workflow {
//     params.sample_list = 'data/samples.dup.tsv'

//     Channel
//         .fromPath(params.sample_list, checkIfExists: true)
//         .set { sample_list }

//     // def sample_list = file('data/samples.dup.tsv')

//     ch_samples_indexed = PREPARE(sample_list).samples

//     ch_samples_consensus = CONSENSUS(ch_samples_indexed).samples

//     ch_samples_consensus.view()

// }
