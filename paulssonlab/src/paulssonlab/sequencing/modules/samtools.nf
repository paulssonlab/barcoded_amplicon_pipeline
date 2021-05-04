include { call_process } from '../functions.nf'

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

def call_SAMTOOLS_SORT(ch) {
    call_process(SAMTOOLS_SORT,
                ch,
                ["bam"],
                [id: { it.bam.name }],
                ["sorted_bam"]) { [it.bam] }
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

def call_SAMTOOLS_INDEX(ch) {
    call_process(SAMTOOLS_INDEX,
                ch,
                ["sorted_bam"],
                [id: { it.sorted_bam.name }],
                ["sorted_bam_index"]) { [it.sorted_bam] }
}
