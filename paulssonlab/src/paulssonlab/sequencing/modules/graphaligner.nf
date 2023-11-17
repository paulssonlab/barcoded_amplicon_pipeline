process GRAPHALIGNER {
    tag "$meta.id"

    time 60.min
    memory 1.GB

    input:
    tuple val(meta), path(reads), path(gfa)

    output:
    tuple val(meta), path("${meta.id}.gaf")

    // TODO
    //conda "${params.conda_env_dir}/graphaligner.yml"

    script:
    """
    GraphAligner -t ${task.cpus} ${meta.graphaligner_args ?: ""} -f ${reads} -g ${gfa} -a ${meta.id}.gaf
    """

    stub:
    """
    touch ${meta.id}.gaf
    """
}
