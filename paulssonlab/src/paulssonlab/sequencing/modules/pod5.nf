process POD5_MERGE {
    tag "$meta.id"
    time { 10.min + 90.min * (task.attempt - 1) }
    memory { 2.5.GB + 6.GB * (task.attempt - 1) }

    input:
    tuple val(meta), path(pod5, stageAs: "pod5/?.pod5")

    output:
    tuple val(meta), path("${meta.id}.pod5")

    conda "${params.conda_env_dir}/pod5.yml"

    script:
    """
    pod5 merge ${meta.pod5_merge_args ?: ""} -o ${meta.id}.pod5 pod5
    """

    stub:
    """
    touch ${meta.id}.pod5
    """
}

process POD5_VIEW {
    tag "$meta.id"
    // time 30.min

    input:
    tuple val(meta), path(pod5, stageAs: "pod5/?.pod5")

    output:
    tuple val(meta), path("view.tsv")

    conda "${params.conda_env_dir}/pod5.yml"

    script:
    """
    pod5 view -t ${task.cpus} ${meta.pod5_view_args ?: ""} -o view.tsv pod5
    """

    stub:
    """
    touch view.tsv
    """
}

process POD5_VIEW_AND_SUBSET {
    tag "$meta.id"
    time 6.hour
    memory { 1.GB + 15.GB * (task.attempt - 1) }

    input:
    tuple val(meta), path(pod5, stageAs: "pod5/?.pod5")

    output:
    tuple val(meta), path("subset/*.pod5")

    conda "${params.conda_env_dir}/pod5.yml"

    script:
    """
    pod5 view -t ${task.cpus} ${meta.pod5_view_args ?: ""} -o view.tsv pod5
    pod5 subset -t ${task.cpus} ${meta.pod5_subset_args ?: ""} -s view.tsv -o subset pod5
    """

    stub:
    """
    mkdir subset
    touch subset/channel-1.pod5
    """
}

process POD5_FILTER {
    tag "$meta.id"
    cpus 2
    time 3.hour
    memory 2.GB

    input:
    tuple val(meta), path(pod5, stageAs: "pod5/?.pod5"), path(read_ids)

    output:
    tuple val(meta), path("${meta.id}.pod5")

    conda "${params.conda_env_dir}/pod5.yml"

    script:
    """
    pod5 filter -t ${task.cpus} ${meta.pod5_filter_args ?: ""} -i ${read_ids} -o ${meta.id}.pod5 pod5
    """

    stub:
    """
    touch ${meta.id}.pod5
    """
}

process SPLIT_READ_IDS {
    tag "$meta.id"
    // time 30.min

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("read_lists/*")

    conda "${params.conda_env_dir}/split_read_ids.yml"

    script:
    """
    ${PIPELINE_ROOT}/paulssonlab/src/paulssonlab/sequencing/bin/split_read_ids.py ${meta.split_read_ids_args ?: ""} ${tsv} read_lists
    """

    stub:
    """
    mkdir read_lists
    touch read_lists/channel-1.tsv
    """
}
