import static functions.*

process CALL_VARIANTS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bam), path(reference)

    output:
    tuple val(meta), path("*.bcf"), path("*.log")

    conda "${params.conda_env_dir}/bcftools.yml"

    script:
    def bcftools_mpileup_args = meta.bcftools_mpileup_args ?: "--max-depth 2000 --max-idepth 2000"
    def bcftools_call_args = meta.bcftools_call_args ?: "--ploidy 1"
    """
    (bcftools mpileup -Ou ${bcftools_mpileup_args} -f ${reference} ${bam} \
     | bcftools call -mv -Ob ${bcftools_call_args} -o variants.bcf) 2> call_variants.log
    """
}

def call_CALL_VARIANTS(ch) {
    call_process(CALL_VARIANTS,
                ch,
                ["bam", "reference", "reference_fai",
                 "bcftools_mpileup_args", "bcftools_call_args"],
                [id: { it.bam.baseName }],
                ["bcf", "call_variants_log"]) { [it.bam, it.reference] }
}

process FILTER_VARIANTS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("*.bcf"), path("*.log")

    conda "${params.conda_env_dir}/bcftools.yml"

    script:
    def bcftools_filter_args = meta.bcftools_filter_args ?: "-i'%QUAL>20'"
    """
    bcftools filter ${bcftools_filter_args} -Ob ${bcf} -o filtered.bcf 2> filter_variants.log
    """
}

def call_FILTER_VARIANTS(ch) {
    call_process(FILTER_VARIANTS,
                ch,
                ["bcf", "bcftools_filter_args"],
                [id: { it.bcf.baseName }],
                ["bcf_filtered", "filter_variants_log"]) { [it.bcf] }
}

process INDEX_VARIANTS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("*.csi")

    conda "${params.conda_env_dir}/bcftools.yml"

    script:
    """
    bcftools index ${bcf}
    """
}

def call_INDEX_VARIANTS(ch) {
    call_process(INDEX_VARIANTS,
                ch,
                ["bcf_filtered"],
                [id: { it.bcf_filtered.baseName }],
                ["bcf_filtered_csi"]) { [it.bcf_filtered] }
}

process GET_CONSENSUS {
    tag "$meta.id"

    input:
    tuple val(meta), path(bcf), path(bcf_csi), path(reference)

    output:
    tuple val(meta), path("consensus.fasta"), path("*.log")

    conda "${params.conda_env_dir}/bcftools.yml"

    script:
    """
    bcftools consensus ${bcf} -f ${reference} -o consensus.fasta 2> get_consensus.log
    """
}

def call_GET_CONSENSUS(ch) {
    call_process(GET_CONSENSUS,
                ch,
                ["bcf_filtered", "bcf_filtered_csi", "reference"],
                [id: { it.bcf_filtered.baseName }],
                ["consensus", "get_consensus_log"]) { [it.bcf_filtered, it.bcf_filtered_csi, it.reference] }
}
