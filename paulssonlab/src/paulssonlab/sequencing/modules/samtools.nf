process SAMTOOLS_FASTQ {
    tag "$meta.id"

    time 10.min
    memory 1.GB

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz")

    conda "${params.conda_env_dir}/samtools.yml"

    script:
    """
    samtools fastq -@ ${task.cpus} ${meta.samtools_fastq_args ?: ""} ${bam} -0 ${meta.id}.fastq.gz
    """

    stub:
    """
    touch ${meta.id}.fastq.gz
    """
}

process SAMTOOLS_IMPORT {
    tag "$meta.id"

    time 10.min
    memory 1.GB

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.bam")

    conda "${params.conda_env_dir}/samtools.yml"

    script:
    """
    samtools import ${meta.samtools_import_args ?: ""} -0 ${fastq} -o ${meta.id}.bam
    """

    stub:
    """
    touch ${meta.id}.bam
    """
}

process SAMTOOLS_MERGE {
    tag "$meta.id"

    time 10.min
    memory 1.GB

    input:
    tuple val(meta), path(input, stageAs: "input/*")

    output:
    tuple val(meta), path("${meta.id}.fastq.gz")

    conda "${params.conda_env_dir}/samtools.yml"

    script:
    """
    samtools merge -@ ${task.cpus} ${meta.samtools_merge_args ?: ""} -o ${meta.id}.fastq.gz ${input}
    """

    stub:
    """
    touch ${meta.id}.fastq.gz
    """
}

process SAMTOOLS_INDEX {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam.name}.bai")

    conda "${params.conda_env_dir}/samtools.yml"

    script:
    """
    samtools index ${bam}
    """

    stub:
    """
    touch ${bam.name}.bai
    """
}

process SAMTOOLS_FAIDX {
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.name}.fai")

    conda "${params.conda_env_dir}/samtools.yml"

    script:
    """
    samtools faidx ${fasta}
    """

    stub:
    """
    touch ${fasta.name}.fai
    """
}
