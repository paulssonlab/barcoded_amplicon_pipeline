import static functions.*

process ANY2FASTA {
    tag "$meta.id"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}.fasta")

    publishDir "${params.output_dir}/${params.references_dir}"

    conda "${params.conda_env_dir}/any2fasta.yml"

    shell:
    '''
    any2fasta -q !{input} | seqkit replace -p '(.*)' -r '!{meta.id}' > !{input.baseName}.fasta
    '''
}

process MERGE_FASTAS {
    tag "$meta.id"

    input:
    tuple val(meta), path("seq")

    output:
    tuple val(meta), path("merged.fasta")

    script:
    """
    cat seq* > merged.fasta
    """
}

process EXTRACT_CONSENSUS {
    // TODO: I think we can remove this (and the below)
    // tag "$meta.id"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("consensus"), path("*.log")

    conda "${params.conda_env_dir}/seqkit.yml"

    script:
    """
    (cat ${fasta} | seqkit replace -p '^(\\S+).*' -r '\$1' \
     | seqkit split - --by-id -f -O consensus) 2> extract_consensus.log
    """
}

def call_EXTRACT_CONSENSUS(ch) {
    call_process(EXTRACT_CONSENSUS,
                ch,
                ["consensus"],
                // [id: { it.consensus.baseName }],
                [:],
                ["consensus_extracted"]) { [it.consensus] }
}
