process SAMTOOLS_SORT {
    tag "$meta.id"

    conda 'envs/mapping.yml'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    shell:
    '''
    samtools sort -@ !{task.cpus} -o !{meta.id}.bam -T !{meta.id} !{bam}
    '''
}

// process SAMTOOLS_INDEX {
//     //
// }
