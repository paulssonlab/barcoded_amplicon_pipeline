import static functions.*

process SAMTOOLS_INDEX {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.name}.bai")

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    """
    samtools index ${bam}
    """
}

def call_SAMTOOLS_INDEX(ch) {
    call_process(SAMTOOLS_INDEX,
                ch,
                ["bam"],
                [id: { it.bam.baseName }],
                ["bam_index"]) { [it.bam] }
}

process SAMTOOLS_FAIDX {
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.name}.fai")

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    """
    samtools faidx ${fasta}
    """
}

def call_SAMTOOLS_FAIDX(ch) {
    call_process(SAMTOOLS_FAIDX,
                ch,
                ["reference"],
                [id: { it.reference.baseName }],
                ["reference_fai"]) { [it.reference] }
}
