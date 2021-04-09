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
    PREPROCESS(samples_in)
    BOWTIE2_BUILD(PREPROCESS.references).index
    PREPROCESS_FINALIZE(PREPROCESS.samples, BOWTIE2_BUILD.index)

    BOWTIE2_INTERLEAVED(samples_in.map { [it, it.reads, it.index] })
        .bam.set { samples }

    emit:
    samples

}

workflow MAIN {
    println SampleSheetParser.load("demo.toml")
    //STAGE_SAMPLE_SHEET()
    //STAGE_REFERENCES(STAGE_READS.samples)
    //ILLUMINA_WHOLE_PLASMID(STAGE_READS.samples)
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
