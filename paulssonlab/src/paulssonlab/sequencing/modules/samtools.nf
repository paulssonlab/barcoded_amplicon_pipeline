import static functions.*

process SAMTOOLS_FASTQ {
    tag "$meta.id"

    time = 10.min
    memory = 1.GB

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.fastq.gz")

    conda "${params.conda_env_dir}/samtools.yml"

    script:
    def args = meta.samtools_fastq_args ?: "-c 1"
    """
    samtools fastq -@ ${task.cpus} ${args} ${bam} -0 ${meta.id}.fastq.gz
    #samtools view -s 17.1 ${bam} | samtools fastq ${args} -0 ${meta.id}.fastq.gz -
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
}

// def call_SAMTOOLS_INDEX(ch) {
//     call_process(SAMTOOLS_INDEX,
//                 ch,
//                 ["bam"],
//                 [id: { it.bam.baseName }],
//                 ["bam_index"]) { [it.bam] }
// }

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
}

// def call_SAMTOOLS_FAIDX(ch) {
//     call_process(SAMTOOLS_FAIDX,
//                 ch,
//                 ["reference"],
//                 [id: { it.reference.baseName }],
//                 ["reference_fai"]) { [it.reference] }
// }
