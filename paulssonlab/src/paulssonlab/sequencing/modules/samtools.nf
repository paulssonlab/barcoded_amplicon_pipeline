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
    #samtools view -s 17.1 ${bam} | samtools fastq ${args} -0 ${meta.id}.fastq.gz -
    """

    stub:
    """
    touch ${meta.id}.fastq.gz
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
    samtools merge -@ ${task.cpus} ${meta.samtools_merge_args ?: ""} ${input} ${meta.id}.fastq.gz
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
