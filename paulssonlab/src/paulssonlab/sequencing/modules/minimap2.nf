import static functions.*

// SEE: https://github.com/nf-core/modules/blob/master/modules/nf-core/minimap2/index/main.nf
process MINIMAP2_INDEX {
    tag "$meta.id"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi")

    conda "${params.conda_env_dir}/minimap2.yml"

    script:
    """
    minimap2 \
        -t ${task.cpus} \\
        -d ${fasta.baseName}.mmi \\
        ${meta.minimap2_index_args ?: ""} \\
        ${fasta}
    """
}

def call_MINIMAP2_INDEX(ch) {
    call_process(MINIMAP2_INDEX,
                 ch,
                 ["reference", "minimap2_index_args"],
                 [id: { it.reference.baseName }],
                 ["index"]) { [it.reference] }
}

// SEE: https://github.com/nf-core/modules/blob/master/modules/nf-core/minimap2/align/main.nf
process MINIMAP2_ALIGN {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("*.bam"), path("*.paf"), path("*.log")

    conda "${params.conda_env_dir}/minimap2.yml"

    script:
    def samtools_command = meta.getOrDefault("sort_bam", true) ? "sort" : "view"
    if (index) {
        reference = index
        def minimap2_align_args = meta.minimap2_align_args ?: "-ax map-ont"
        def output = "| samtools ${samtools_command} --threads ${task.cpus} -o ${meta.id}.bam -"
    }
    else {
        reference = reads
        def minimap2_align_args = meta.minimap2_align_args ?: "-x ava-ont"
        def output = "-o ${meta.id}.paf"
    }
    """
    (minimap2 \\
        ${meta.minimap2_align_args ?: ""} \\
        -t ${task.cpus} \\
        "${reference}" \\
        "${reads}" \\
        ${output}) 2> ${meta.id}.minimap2.log
    """
}

def call_MINIMAP2_ALIGN(ch) {
    call_process(MINIMAP2_ALIGN,
                 ch,
                 ["reads", "index", "minimap2_align_args", "sort_bam"],
                 [id: { it.reads.baseName }],
                 ["bam", "paf", "minimap2_log"]) { [it.reads, it.index] }
}
