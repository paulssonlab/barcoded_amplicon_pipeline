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

    stub:
    """
    touch ${model}
    """
}

process DORADO_BASECALLER {
    tag "$meta.id"
    label "dorado_gpu"
    cpus 2
    memory 20.GB
    time 90.min // TODO: adjust based on total input file size

    input:
    tuple val(meta), path("pod5/?.pod5"), path(dorado_model)

    output:
    tuple val(meta), path("${meta.id}.bam")

    // TODO
    // conda "${params.conda_env_dir}/dorado.yml"

    script:
    """
    module load gcc/9.2.0
    module load cuda/12.1
    # TODO
    echo -n "HOSTNAME: "
    hostname -a
    echo -n "GPU: "
    nvidia-smi --query-gpu=name --format=csv,noheader
    dorado basecaller ${meta.dorado_basecaller_args ?: ""} ${dorado_model} pod5 > ${meta.id}.bam
    """

    stub:
    """
    touch ${meta.id}.bam
    """
}

process DORADO_DUPLEX {
    tag "$meta.id"
    label "dorado_gpu"
    cpus 2
    memory 20.GB
    time 90.min // TODO: adjust based on total input file size

    input:
    tuple val(meta), path("pod5/?.pod5"), path(dorado_model), path(dorado_duplex_model)

    output:
    tuple val(meta), path("${meta.id}.bam")

    // TODO
    // conda "${params.conda_env_dir}/dorado.yml"

    script:
    """
    module load gcc/9.2.0
    module load cuda/12.1
    # TODO
    echo -n "HOSTNAME: "
    hostname -a
    echo -n "GPU: "
    nvidia-smi --query-gpu=name --format=csv,noheader
    dorado duplex ${meta.dorado_duplex_args ?: ""} ${dorado_model} pod5 > ${meta.id}.bam
    """

    stub:
    """
    touch ${meta.id}.bam
    """
}

process DORADO_DUPLEX_WITH_PAIRS {
    tag "$meta.id"
    label "dorado_duplex_only_gpu"
    cpus 2
    memory 16.GB
    time 90.min // TODO: adjust based on total input file size

    input:
    tuple val(meta), path("pod5/?.pod5"), path(pairs), path(dorado_model), path(dorado_duplex_model)

    output:
    tuple val(meta), path("${meta.id}_duplex.bam")

    // TODO
    // conda "${params.conda_env_dir}/dorado.yml"

    script:
    """
    module load gcc/9.2.0
    module load cuda/12.1
    # TODO
    echo -n "HOSTNAME: "
    hostname -a
    echo -n "GPU: "
    nvidia-smi --query-gpu=name --format=csv,noheader
    dorado duplex --pairs ${pairs} ${meta.dorado_duplex_args ?: ""} ${dorado_model} pod5 > ${meta.id}_duplex.bam
    """

    stub:
    """
    touch ${meta.id}_duplex.bam
    """
}
