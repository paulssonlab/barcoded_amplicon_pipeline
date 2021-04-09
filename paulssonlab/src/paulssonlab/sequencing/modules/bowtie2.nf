process BOWTIE2_BUILD {
    tag "$meta"

    conda 'envs/mapping.yml'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2'), emit: index

    script:
    '''
    mkdir bowtie2
    bowtie2-build --threads !{task.cpus} !{fasta} bowtie2/!{fasta.baseName}
    '''
}

process BOWTIE2_INTERLEAVED {
    tag "$meta.id"

    conda 'envs/mapping.yml'

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.log'), emit: log

    shell:
    '''
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
    (bowtie2 \
        --threads !{task.cpus} \
        -x $INDEX --interleaved !{reads} \
        | samtools view -@ !{task.cpus} -Sbh -o !{meta.id}.bam -) 2> !{meta.id}.bowtie.log
    '''
}
