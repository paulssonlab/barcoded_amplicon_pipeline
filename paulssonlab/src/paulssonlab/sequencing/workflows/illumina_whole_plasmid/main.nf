include { PREPARE_SAMPLE_SHEET;
          PREPARE_READS;
          PREPARE_REFERENCES;
          MERGE_INDICES } from '../prepare.nf'

include { BOWTIE2_BUILD;
          BOWTIE2_INTERLEAVED } from '../../modules/bowtie2.nf'

workflow ILLUMINA_WHOLE_PLASMID {
    take:
    samples_in

    main:
    PREPARE_REFERENCES(samples_in)
    // PREPARE_REFERENCES.references.view()
    // PREPARE_REFERENCES.samples.view()
    // BOWTIE2_BUILD(PREPARE_REFERENCES.references).index
    // MERGE_INDICES(PREPARE_REFERENCES.samples, BOWTIE2_BUILD.index)

    // BOWTIE2_INTERLEAVED(samples_in.map { [it, it.reads, it.index] })
    //     .bam.set { samples }

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
