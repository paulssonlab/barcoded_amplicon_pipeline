process JOIN_GAF {
    tag "$meta.id"



    input:
    tuple val(meta), path(input), path(gaf)

    output:
    tuple val(meta), path("${meta.id}.${meta.tabular_format}")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/join_gaf.yml"

    script:
    def input_format = meta.join_gaf_input_format ?: "arrow"
    def output_format = meta.join_gaf_output_format ?: "arrow"
    """
    bin/join_gaf.py ${meta.join_gaf_args ?: ""} \
        --input-format ${meta.tabular_format} \
        --output-format ${meta.tabular_format} \
        --gaf ${gaf} \
        ${input} \
        ${meta.id}.${meta.tabular_format}
    """

    stub:
    """
    touch ${meta.id}.${meta.tabular_format}
    """
}
