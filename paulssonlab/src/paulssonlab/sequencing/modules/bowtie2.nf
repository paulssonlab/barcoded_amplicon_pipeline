import static functions.*

// SEE: https://github.com/nf-core/modules/blob/master/modules/bowtie2/build/main.nf
process BOWTIE2_BUILD {
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2')

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    """
    mkdir bowtie2
    bowtie2-build --threads ${task.cpus} ${meta.bowtie2_build_args ?: ""} ${fasta} bowtie2/${fasta.baseName}
    """
}

def call_BOWTIE2_BUILD(ch) {
    call_process(BOWTIE2_BUILD,
                 ch,
                 ["reference", "bowtie2_build_args"],
                 [id: { it.reference.baseName }],
                 ["index"]) { [it.reference] }
}

// SEE: https://github.com/nf-core/modules/blob/master/modules/bowtie2/align/main.nf
process BOWTIE2_ALIGN {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('*.bam'), path('*.log')

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    // TODO: copy nf-core to implement non-interleaved modes: -1 aaa_1.fastq -2 aaa_2.fastq / -U aaa.fastq
    def bowtie2_args = meta.bowtie2_args ?: "--end-to-end"
    def samtools_command = meta.getOrDefault("sort_bam", true) ? "sort" : "view"
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/.rev.1.bt2//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/.rev.1.bt2l//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1

    (bowtie2 \
        --threads ${task.cpus} \
        ${bowtie2_args} \
        -x \$INDEX --interleaved ${reads} \
        | samtools ${samtools_command} --threads ${task.cpus} -o ${meta.id}.bam -) 2> ${meta.id}.bowtie2.log
    """
}

def call_BOWTIE2_ALIGN(ch) {
    call_process(BOWTIE2_ALIGN,
                 ch,
                 ["reads", "index", "bowtie2_args", "sort_bam"],
                 [id: { it.reads.baseName }],
                 ["bam", "bowtie2_log"]) { [it.reads, it.index] }
}
