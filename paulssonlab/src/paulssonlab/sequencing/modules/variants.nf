import static functions.*

process CALL_VARIANTS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(reference)

    output:
    tuple val(meta), path("*.bcf"), path("*.log")

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    def bcftools_mpileup_args = meta.bcftools_mpileup_args ?: "--max-depth 2000 --max-idepth 2000"
    def bcftools_call_args = meta.bcftools_call_args ?: "--ploidy 1"
    """
    (bcftools mpileup -Ou ${bcftools_mpileup_args} -f ${reference} ${bam}
     | bcftools call -mv -Ob ${bcftools_call_args} -o ${meta.id}.bcf) 2> ${meta.id}.call_variants.log
    """
}

process FILTER_VARIANTS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("*.bcf"), path("*.log")

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    def bcftools_filter_args = meta.bcftools_filter_args ?: "-i'%QUAL>20'"
    """
    bcftools filter ${bcftools_filter_args} -Ob ${bcf} -o ${meta.id}.filtered.bcf 2> ${meta.id}.filter_variants.log
    """
}

// process INDEX_VARIANTS {
//     tag "$meta.id"

//     input:
//     tuple val(meta), path(bcf)

//     output:
//     tuple val(meta), path("*.csi")

//     conda "${params.conda_env_dir}/mapping.yml"

//     shell:
//     '''
//     bcftools index ${bcf}
//     '''
// }

process GET_CONSENSUS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bcf), path(csi), path(reference)

    output:
    tuple val(meta), path("*.fasta"), path("*.log")

    conda "${params.conda_env_dir}/mapping.yml"

    script:
    """
    bcftools consensus ${bcf} -f ${reference} -o ${meta.id}.consensus.fasta 2> ${meta.id}.get_consensus.log
    """
}
