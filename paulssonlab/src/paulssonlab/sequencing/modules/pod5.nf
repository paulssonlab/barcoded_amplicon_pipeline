process POD5_MERGE {
    tag "$meta.id"
    label "pod5"

    input:
    tuple val(meta), path(pod5, stageAs: "?.pod5")

    output:
    tuple val(meta), path("merged.pod5")

    conda "${params.conda_env_dir}/pod5.yml"

    script:
    """
    pod5 merge ${meta.pod5_merge_args ?: ""} -o merged.pod5 ${pod5}
    """
}

process POD5_VIEW {
    tag "$meta.id"
    label "pod5"

    input:
    tuple val(meta), path(pod5, stageAs: "?.pod5")

    output:
    tuple val(meta), path("view.tsv")

    conda "${params.conda_env_dir}/pod5.yml"

    script:
    """
    pod5 view -t ${task.cpus} ${meta.pod5_view_args ?: ""} -o view.tsv ${pod5}
    """
}

process POD5_FILTER {
    tag "$meta.id"
    label "pod5"

    input:
    tuple val(meta), path(pod5, stageAs: "?.pod5"), path(read_ids)

    output:
    tuple val(meta), path("${meta.id}.pod5")

    conda "${params.conda_env_dir}/pod5.yml"

    script:
    """
    pod5 filter -t ${task.cpus} ${meta.pod5_filter_args ?: ""} -i ${read_ids} -o ${meta.id}.pod5 ${pod5}
    """
}

process SPLIT_READ_IDS {
    tag "$meta.id"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("read_lists/*")

    conda "${params.conda_env_dir}/split_read_ids.yml"

    script:
    """
    ${System.env['src']}/sequencing/bin/split_read_ids.py ${meta.split_read_ids_args ?: ""} ${tsv} read_lists
    """
}
