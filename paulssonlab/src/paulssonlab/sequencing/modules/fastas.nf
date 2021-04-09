process ANY2FASTA {
    tag "$meta.id"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}.fasta")

    publishDir 'test/references'

    conda 'envs/any2fasta.yml'

    shell:
    '''
    any2fasta -q !{input} | seqkit replace -p '(.*)' -r '!{meta.id}' > !{meta.id}.fasta
    '''
}

process MERGE_FASTAS {
    tag "$meta"

    input:
    tuple val(meta), path('seq')

    output:
    tuple val(meta), path('reference.fasta')

    shell:
    '''
    cat seq* > reference.fasta
    '''
}

// process EXTRACT_CONSENSUS {
//     //
// }
