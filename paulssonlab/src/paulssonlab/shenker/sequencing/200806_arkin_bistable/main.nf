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

// include { GET_REGISTRY_SEQS;
//           ANY2FASTA;
//           MERGE_FASTAS } from './modules.nf'

def zip(ch_a, ch_b) {
    //return ch_a.map { [it] }.combine(ch_b.map { [it] }).map { it.transpose().collectEntries() }
    return ch_a.map { [it] }.combine(ch_b.map { [it] }).map { it.transpose() }
}

process GET_REGISTRY_SEQS {
    input:
    stdin

    output:
    path('*.gb', type: 'file')

    storeDir 'test/references'

    shell:
    '''
    #!/usr/bin/env python
    import sys
    import shutil
    import re

    SOURCE = "/Users/jacob/Dropbox (Personal)/Research/Paulsson/paulssonlab/paulssonlab/src/paulssonlab/shenker/sequencing/200806_arkin_bistable/references/pLIB219.gb"

    ids = re.split(r"\\s*,\\s*", sys.stdin.read().rstrip())
    for i, id in enumerate(ids):
        shutil.copy(SOURCE, f"{id}.gb")
        #with open(f"{id}.gb", "w") as f:
        #    f.write(f"BLAH {id} BLEE\\n")
    '''
}

process ANY2FASTA {
    tag "$meta.id"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}.fasta")

    publishDir 'test/references'

    conda 'envs/any2fasta.yml'

    shell:
    """
    any2fasta -q $input | seqkit replace -p '(.*)' -r '${meta.id}' > ${meta.id}.fasta
    """
}

process MERGE_FASTAS {
    input:
    tuple val(meta), path('seq')

    output:
    path('reference.fasta')

    shell:
    """
    cat seq* > reference.fasta
    """
}

workflow {
    params.samples_list = 'data/samples.dup.tsv'

    Channel
        .fromPath(params.samples_list, checkIfExists: true)
        .splitCsv(sep:'\t')
        .map { row -> [id: row[0], references: row[1].split('\s*,\s*')] }
        .set { samples }

    samples
        .flatMap { meta -> meta.references }
        .unique()
        .collect()
        .set { reference_names }

    // zip(reference_names, GET_REGISTRY_SEQS(reference_names.map { it.join(",") }))
    //     .set { references }

    GET_REGISTRY_SEQS(reference_names.map { it.join(",") })
        //.map { it.collect { f -> [[id: f.getBaseName()], f] }}
        .flatMap { it.collect { f -> [[id: f.getBaseName()], f] }}
        .set { ch_references_orig_format }

    ANY2FASTA(ch_references_orig_format).view()

    //ANY2FASTA().view()

    //references.view()

    //samples.view()

    // samples
    //     .map { [it, file("data/${it[0]}\.fastq")] }.view()
    //     .set { reads }

}
