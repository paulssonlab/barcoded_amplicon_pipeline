process SAMTOOLS_SORT {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.baseName}.sorted.bam")

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    """
    samtools sort -@ ${task.cpus} -o ${bam.baseName}.sorted.bam ${bam}
    """
}

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
