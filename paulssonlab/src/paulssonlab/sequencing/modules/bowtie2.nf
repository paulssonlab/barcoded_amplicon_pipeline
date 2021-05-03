// SEE: https://github.com/nf-core/modules/blob/master/software/bowtie2/build/main.nf
process BOWTIE2_BUILD {
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2'), emit: index

    conda "${params.conda_env_dir}/mapping.yml"

    shell:
    '''
    mkdir bowtie2
    bowtie2-build --threads !{task.cpus} !{fasta} bowtie2/!{fasta.baseName}
    '''
}

// SEE: https://github.com/nf-core/modules/blob/master/software/bowtie2/align/main.nf
process BOWTIE2_INTERLEAVED {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('*.bam'), path('*.log')

    conda "${params.conda_env_dir}/mapping.yml"

    shell:
    '''
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
    (bowtie2 \
        --threads !{task.cpus} \
        -x $INDEX --interleaved !{reads} \
        | samtools view -@ !{task.cpus} -Sbh -o !{meta.id}.bam -) 2> !{meta.id}.bowtie.log
    '''
}
