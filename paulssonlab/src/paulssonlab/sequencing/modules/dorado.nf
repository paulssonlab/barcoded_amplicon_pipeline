process DORADO_DOWNLOAD {
    tag "$meta.id"
    label "local"

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
    label "dorado_gpu"
    cpus = 2 // TODO: useful?
    time = 12.hours // TODO: adjust based on total input file size
    memory = 32.GB
    scratch true
    stageInMode "copy"

    input:
    tuple val(meta), path("pod5/?.pod5"), path(dorado_model), path(dorado_duplex_model)

    output:
    tuple val(meta), path("${meta.id}.bam")

    // TODO
    // conda "${params.conda_env_dir}/dorado.yml"

    script:
    """
    # TODO
    echo -n HOSTNAME:
    hostname -a
    module load gcc/9.2.0
    module load cuda/11.7
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
    # TODO
    module load gcc/9.2.0
    module load cuda/11.7
    dorado basecaller ${meta.dorado_basecaller_args ?: ""} ${dorado_model} pod5 > ${meta.id}.bam
    """
}
