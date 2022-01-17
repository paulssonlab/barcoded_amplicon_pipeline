include { call_process } from '../../functions.nf'

include { PREPARE_SAMPLE_SHEET;
          PREPARE_READS;
          PREPARE_REFERENCES;
          MERGE_INDEXES } from '../prepare.nf'

include { BOWTIE2_BUILD;
          call_BOWTIE2 } from '../../modules/bowtie2.nf'

include { call_SAMTOOLS_SORT;
          call_SAMTOOLS_INDEX } from '../../modules/samtools.nf'

workflow ILLUMINA_WHOLE_PLASMID {
    take:
    samples_in

    main:
    PREPARE_REFERENCES(samples_in)

    // [reads_path:4309_APA4309_321712w_AE5.fastq, param_set:default, reads_prefix:4309_APA4309_321712w_AE5, reference_names:[pLIB219, pLIB220], name:4309_APA4309_321712w_AE5, id:0, run_path:default/4309_APA4309_321712w_AE5, reads:/tmp/paulssonlab-sequencing/200806_arkin_bistable_sequencing/data/4309_APA4309_321712w_AE5.fastq, references:[/tmp/paulssonlab-sequencing/200806_arkin_bistable_sequencing/data/references/pLIB219.gb, /tmp/paulssonlab-sequencing/200806_arkin_bistable_sequencing/data/references/pLIB220.gb]]
    // [reads_path:4310_APA4310_321712w_AE6.fastq, param_set:default, reads_prefix:4310_APA4310_321712w_AE6, reference_names:[pLIB219, pLIB222], name:4310_APA4310_321712w_AE6, id:1, run_path:default/4310_APA4310_321712w_AE6, reads:/tmp/paulssonlab-sequencing/200806_arkin_bistable_sequencing/data/4310_APA4310_321712w_AE6.fastq, references:[/tmp/paulssonlab-sequencing/200806_arkin_bistable_sequencing/data/references/pLIB219.gb, /tmp/paulssonlab-sequencing/200806_arkin_bistable_sequencing/data/references/pLIB222.gb]]
    // [reads_path:4311_APA4311_321712w_AE7.fastq, param_set:default, reads_prefix:4311_APA4311_321712w_AE7, reference_names:[pLIB221, pLIB222], name:4311_APA4311_321712w_AE7, id:2, run_path:default/4311_APA4311_321712w_AE7, reads:/tmp/paulssonlab-sequencing/200806_arkin_bistable_sequencing/data/4311_APA4311_321712w_AE7.fastq, references:[/tmp/paulssonlab-sequencing/200806_arkin_bistable_sequencing/data/references/pLIB221.gb, /tmp/paulssonlab-sequencing/200806_arkin_bistable_sequencing/data/references/pLIB222.gb]]

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
