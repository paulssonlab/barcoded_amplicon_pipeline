#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

include { GET_REGISTRY_SEQS;
          ANY2FASTA;
          MERGE_FASTAS;
          BOWTIE2_BUILD } from './modules.nf'

workflow {
    params.samples_list = 'data/samples.dup.tsv'

    Channel
        .fromPath(params.samples_list, checkIfExists: true)
        .splitCsv(sep:'\t')
        .map { row -> [id: row[0], references: row[1].split('\s*,\s*') as Set] }
        .map { [reads: file("data/${it.id}.fastq"), *:it]}
        .set { ch_samples }

    ch_samples
        .flatMap { meta -> meta.references }
        .unique()
        .collect()
        .set { ch_reference_names }

    GET_REGISTRY_SEQS(ch_reference_names.map { it.join(",") })
        .flatMap { it.collect { f -> [[id: f.getBaseName()], f] }}
        .set { ch_references_orig_format }

    ANY2FASTA(ch_references_orig_format)
        .set { ch_references }

    ch_references
        .map { [[it[0].id, it[1]]] }
        .collect()
        .map { it.collectEntries() }
        .set { ch_collected_references }

    ch_samples
        .map { it.references }
        .unique()
        .set { ch_reference_sets }

    ch_reference_sets
        .combine(ch_collected_references)
        .map { ref_set, refs -> [ref_set, ref_set.collect { refs[it] }] }
        .set { ch_reference_set_fastas }

    MERGE_FASTAS(ch_reference_set_fastas)
        .set { ch_merged_references }

    BOWTIE2_BUILD(ch_merged_references).index
        .set { ch_indexes }

    ch_samples
        .map { [it.references, it] }
        .join(ch_indexes)
        .map { [index: it[2], *:it[1]] }
        .view()

    //ANY2FASTA().view()

    //references.view()

    //samples.view()

    // samples
    //     .map { [it, file("data/${it[0]}\.fastq")] }.view()
    //     .set { reads }

}
