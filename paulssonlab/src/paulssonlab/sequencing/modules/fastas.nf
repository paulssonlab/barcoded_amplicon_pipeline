import java.nio.file.Paths

process ANY2FASTA {
    tag "$meta.id"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}.fasta")

    publishDir "${workDir}/${params.references_dir}"

    conda "${params.conda_env_dir}/any2fasta.yml"

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
