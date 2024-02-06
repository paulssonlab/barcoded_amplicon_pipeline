process FIND_DUPLEX_PAIRS {
    tag "$meta.id"

    time 10.min
    memory 4.GB

    input:
    tuple val(meta), path(gfa), path(bam), path(gaf)

    output:
    tuple val(meta), path("${meta.id}_pairs.txt")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/scripts.yml"

    script:
    """
    ${src}/sequencing/bin/find_duplex_pairs.py ${meta.find_duplex_reads_args ?: ""} --gfa ${gfa} --gaf ${gaf} ${bam} ${meta.id}_pairs.txt
    """

    stub:
    """
    touch ${meta.id}_pairs.txt
    """
}

process JOIN_GAF {
    tag "$meta.id"

    time { 10.min + 40.min * (task.attempt - 1) }
    memory { 4.GB + 4.GB * (task.attempt - 1) }

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

    time { 10.min + 50.min * (task.attempt - 1) }
    memory { 8.GB + 16.GB * (task.attempt - 1) }

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

    time 180.min
    memory { 8.GB + 16.GB * (task.attempt - 1) }

    input:
    tuple val(meta), path(input, stageAs: "input/*")

    output:
    tuple val(meta), path("${meta.id}.${meta.output_format}")

    // TODO
    conda "/home/jqs1/micromamba/envs/medaka"
    // conda "${params.conda_env_dir}/scripts.yml"

    script:
    """
    ${src}/sequencing/bin/consensus.py --skip-consensus ${meta.prepare_consensus_args ?: ""} --group ${meta.group} --input-format ${meta.input_format} --output-format ${meta.output_format} --output ${meta.id}.${meta.output_format} ${input}
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

    time 180.min
    memory { 8.GB + 16.GB * (task.attempt - 1) }

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

     time 180.min
    memory { 8.GB + 16.GB * (task.attempt - 1) }

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

    time { 120.min + 120.min * (task.attempt - 1) }
    memory { 1.GB + 7.GB * (task.attempt - 1) }

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

    time { 10.min + 60.min * (task.attempt - 1) }
    memory { 2.GB + 6.GB * (task.attempt - 1) }

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
