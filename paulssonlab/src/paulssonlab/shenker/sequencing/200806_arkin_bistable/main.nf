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
    tuple val(meta), path('reference.fasta')

    shell:
    """
    cat seq* > reference.fasta
    """
}

// process BOWTIE2_BUILD {
//     conda: '/envs/mapping.yml'
// }

process BOWTIE2_BUILD {
    tag "$meta"

    conda 'envs/mapping.yml'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2'), emit: index

    script:
    """
    mkdir bowtie2
    bowtie2-build --threads $task.cpus $fasta bowtie2/${fasta.baseName}
    """
}

// process BOWTIE2_INTERLEAVED {
//     //
// }

workflow {
    params.samples_list = 'data/samples.dup.tsv'

    Channel
        .fromPath(params.samples_list, checkIfExists: true)
        .splitCsv(sep:'\t')
        .map { row -> [id: row[0], references: row[1].split('\s*,\s*') as Set] }
        //.map { row -> [id: row[0], references: row[1].split('\s*,\s*')] }
        .set { ch_samples }

    ch_samples
        .flatMap { meta -> meta.references }
        .unique()
        .collect()
        .set { ch_reference_names }

    // zip(reference_names, GET_REGISTRY_SEQS(reference_names.map { it.join(",") }))
    //     .set { references }

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
        //.map { it.references as Set }
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

    //ch_indexes.view()

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
