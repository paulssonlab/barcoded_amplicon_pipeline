process GRAPHALIGNER {
    tag "$meta.id"

    time 60.min
    memory 1.GB

    input:
    tuple val(meta), path(reads), path(gfa)

    output:
    tuple val(meta), path("${meta.id}.gaf")

    conda "${params.conda_env_dir}/graphaligner.yml"

    script:
    if (meta.getOrDefault("graphaligner_strip_tags", true)) {
        """
        seqkit replace -p "\\s.+" -o ${meta.id}_notags.fastq.gz ${reads}
        GraphAligner -t ${task.cpus} ${meta.graphaligner_args ?: ""} -f ${meta.id}_notags.fastq.gz -g ${gfa} -a ${meta.id}.gaf
        """
    } else {
        """
        GraphAligner -t ${task.cpus} ${meta.graphaligner_args ?: ""} -f ${reads} -g ${gfa} -a ${meta.id}.gaf
        """
    }

    stub:
    """
    touch ${meta.id}.gaf
    """
}
