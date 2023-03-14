import static functions.*

process ANY2FASTA {
    tag "$meta.id"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}.fasta")

    publishDir params.references_dir, mode: "copy"

    conda "${params.conda_env_dir}/any2fasta.yml"

    script:
    """
    any2fasta -q ${input} | seqkit replace -p '(.*)' -r '${meta.id}' > ${input.baseName}.fasta
    """
}

process MERGE_FILES {
    tag "$meta.id"

    input:
    tuple val(meta), path("input"), val(ext)

    output:
    tuple val(meta), path("merged.${ext}")

    // SEE: https://github.com/nextflow-io/nextflow/discussions/2813
    // MERGE_FILES is run for each unique reference list, not for each
    // run_path, so we need a separate process to publish merged fastas
    // publishDir { meta.run_output_dir }, mode: "copy"

    script:
    """
    cat input* > merged.${ext}
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
                [id: { it.consensus.baseName }],
                ["consensus_extracted"]) { [it.consensus] }
}
