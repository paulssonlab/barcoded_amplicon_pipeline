include { call_process } from '../functions.nf'

// SEE: https://github.com/nf-core/modules/blob/master/software/bowtie2/build/main.nf
process BOWTIE2_BUILD {
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2'), emit: index

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    def bowtie2_build_args = meta.bowtie2_build_args ?: ""
    """
    mkdir bowtie2
    bowtie2-build --threads ${task.cpus} ${bowtie2_build_args} ${fasta} bowtie2/${fasta.baseName}
    """
}

// SEE: https://github.com/nf-core/modules/blob/master/software/bowtie2/align/main.nf
process BOWTIE2 {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('*.bam'), path('*.log')

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    // TODO: implement non-interleaved modes: -1 aaa_1.fastq -2 aaa_2.fastq / -U aaa.fastq
    def bowtie2_args = meta.bowtie2_args ?: "--end-to-end"
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
    (bowtie2 \
        --threads ${task.cpus} \
        ${bowtie2_args} \
        -x \$INDEX --interleaved ${reads} \
        | samtools view -@ ${task.cpus} -Sbh -o ${meta.id}.bam -) 2> ${meta.id}.bowtie2.log
    """
}

def call_BOWTIE2(ch) {
    call_process(BOWTIE2,
                 ch,
                 ["reads", "index", "bowtie2_args"],
                 [id: { it.reads.baseName }],
                 ["bam", "bowtie2_log"]) { [it.reads, it.index] }
}
