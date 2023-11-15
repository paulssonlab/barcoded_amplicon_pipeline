process JOIN_GAF {
    tag "$meta.id"

    time 10.min
    memory 8.GB

    input:
    tuple val(meta), path(input, stageAs: "input/*"), path(gaf)

    output:
    tuple val(meta), path("${meta.id}.${meta.output_format}")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/scripts.yml"

    script:
    """
    ${src}/sequencing/bin/join_gaf.py ${meta.join_gaf_args ?: ""} --input-format ${meta.input_format} --output-format ${meta.output_format} --gaf ${gaf} ${input} ${meta.id}.${meta.output_format}
    """

    stub:
    """
    touch ${meta.id}.${meta.output_format}
    """
}

process PREPARE_READS {
    tag "$meta.id"

    time 10.min
    memory 12.GB

    input:
    tuple val(meta), path(input, stageAs: "input/*"), path(gfa)

    output:
    tuple val(meta), path("${meta.id}.${meta.output_format}")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/scripts.yml"

    script:
    """
    ${src}/sequencing/bin/prepare_reads.py ${meta.prepare_reads_args ?: ""} --input-format ${meta.input_format} --output-format ${meta.output_format} --gfa ${gfa} ${input} ${meta.id}.${meta.output_format}
    """

    stub:
    """
    touch ${meta.id}.${meta.output_format}
    """
}

process PREPARE_CONSENSUS {
    tag "$meta.id"

    time 4.hours
    memory 8.GB

    errorStrategy "retry"

    input:
    tuple val(meta), path(input, stageAs: "input/*")

    output:
    tuple val(meta), path("${meta.id}.${meta.output_format}")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/scripts.yml"

    script:
    """
    ${src}/sequencing/bin/consensus.py ${meta.consensus_prepare_args ?: ""} --group ${meta.group} --input-format ${meta.input_format} --output-format ${meta.output_format} --output ${meta.id}.${meta.output_format} ${input}
    """

    stub:
    """
    touch ${meta.id}.${meta.output_format}
    """
}

// TODO: consolidate with CONSENSUS (maybe CONSENSUS_PREPARED?)
// so that one process definition can work for all two/three uses
process CONSENSUS_PREPARED {
    tag "$meta.id"

    time 2.hours
    memory 4.GB

    errorStrategy "retry"

    input:
    tuple val(meta), path(input, stageAs: "input/*")

    output:
    tuple val(meta), path("${meta.id}.${meta.output_format}"), path("${meta.id}.fasta")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/scripts.yml"

    script:
    """
    ${src}/sequencing/bin/consensus.py ${meta.consensus_args ?: ""} --input-format ${meta.input_format} --output-format ${meta.output_format} --output ${meta.id}.${meta.output_format} --fasta ${meta.id}.fasta ${input}
    """

    stub:
    """
    touch ${meta.id}.${meta.output_format}
    touch ${meta.id}.fasta
    """
}

process CONSENSUS {
    tag "$meta.id"

    time 4.hours
    memory 8.GB

    errorStrategy "retry"

    input:
    tuple val(meta), path(input, stageAs: "input/*")

    output:
    tuple val(meta), path("${meta.id}.${meta.output_format}"), path("${meta.id}.fasta")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/scripts.yml"

    script:
    """
    ${src}/sequencing/bin/consensus.py ${meta.consensus_args ?: ""} --group ${meta.group} --input-format ${meta.input_format} --output-format ${meta.output_format} --output ${meta.id}.${meta.output_format} --fasta ${meta.id}.fasta ${input}
    """

    stub:
    """
    touch ${meta.id}.${meta.output_format}
    touch ${meta.id}.fasta
    """
}

process REALIGN {
    tag "$meta.id"

    time 60.min
    memory 1.GB

    input:
    tuple val(meta), path(input, stageAs: "input/*"), path(gfa)

    output:
    tuple val(meta), path("${meta.id}.${meta.output_format}")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/scripts.yml"

    script:
    """
    ${src}/sequencing/bin/realign.py ${meta.realign_args ?: ""} --input-format ${meta.input_format} --output-format ${meta.output_format} --gfa ${gfa} ${input} ${meta.id}.${meta.output_format}
    """

    stub:
    """
    touch ${meta.id}.${meta.output_format}
    """
}

process EXTRACT_SEGMENTS {
    tag "$meta.id"

    time 10.min
    memory 1.GB

    input:
    tuple val(meta), path(input, stageAs: "input/*"), path(gfa)

    output:
    tuple val(meta), path("${meta.id}.${meta.output_format}")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/scripts.yml"

    script:
    """
    ${src}/sequencing/bin/extract_segments.py ${meta.extract_segments_args ?: ""} --input-format ${meta.input_format} --output-format ${meta.output_format} --gfa ${gfa} ${input} ${meta.id}.${meta.output_format}
    """

    stub:
    """
    touch ${meta.id}.${meta.output_format}
    """
}
