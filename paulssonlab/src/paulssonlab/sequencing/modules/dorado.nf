process DORADO_DOWNLOAD {
    tag "$meta.id"

    input:
    tuple val(meta), val(model)

    output:
    tuple val(meta), path(model)

    // TODO
    // conda "${params.conda_env_dir}/dorado.yml"

    script:
    """
    dorado download ${meta.dorado_download_args ?: ""} --model ${model}
    """
}

process DORADO_DUPLEX {
    tag "$meta.id"

    input:
    tuple val(meta), path("pod5/?.pod5"), path(dorado_model), path(dorado_duplex_model)

    output:
    tuple val(meta), path("${meta.id}.bam")

    // TODO
    // conda "${params.conda_env_dir}/dorado.yml"

    script:
    """
    dorado duplex ${meta.dorado_duplex_args ?: ""} ${dorado_model} pod5 > ${meta.id}.bam
    """
}

process DORADO_BASECALLER {
    tag "$meta.id"

    input:
    tuple val(meta), path("pod5/?.pod5"), path(dorado_model)

    output:
    tuple val(meta), path("${meta.id}.bam")

    // TODO
    // conda "${params.conda_env_dir}/dorado.yml"

    script:
    """
    dorado basecaller ${meta.dorado_basecaller_args ?: ""} ${dorado_model} pod5 > ${meta.id}.bam
    """
}
