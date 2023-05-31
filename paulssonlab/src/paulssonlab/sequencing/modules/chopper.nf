import static functions.*

process CHOPPER {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("filtered.fastq.gz")

    conda "${params.conda_env_dir}/chopper.yml"

    script:
    def chopper_args = []
    if (meta.min_read_length) {
        chopper_args << "--minlength ${meta.min_read_length}"
    }
    if (meta.max_read_length) {
        chopper_args << "--maxlength ${meta.max_read_length}"
    }
    if (meta.min_read_quality) {
        chopper_args << "--quality ${meta.min_read_quality}"
    }
    if (meta.read_head_crop) {
        chopper_args << "--headcrop ${meta.read_head_crop}"
    }
    if (meta.read_tail_crop) {
        chopper_args << "--tailcrop ${meta.read_tail_crop}"
    }
    if (meta.decontaminate_reference) {
        chopper_args << "--contam ${meta.decontaminate_reference}"
    }
    chopper_args = chopper_args.join(" ")
    """
    gunzip -c ${reads} | chopper --threads ${task.cpus} ${chopper_args} | gzip > filtered.fastq.gz
    """
}

def call_CHOPPER(ch) {
    def chopper_params = ["min_read_length", "max_read_length", "min_read_quality",
                          "read_head_crop", "read_tail_crop", "decontaminate_reference"]
    call_process(CHOPPER,
                 ch,
                 ["reads", *chopper_params],
                 [id: { it.reads?[0]?.baseName }],
                 ["filtered_reads"], { chopper_params.collect { p -> it[p] }.any() }) { [it.reads] }
        .map { rename_key(it, "filtered_reads", "filtered_reads") { reads -> [reads] } }
}
