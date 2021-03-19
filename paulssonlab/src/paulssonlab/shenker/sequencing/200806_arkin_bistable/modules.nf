process SCP {
    input:
    val(remote_path)

    output:
    path('samples.tsv')

    shell:
    '''
    scp !{remote_path} !{}
    '''
}

process BOWTIE2_BUILD {
    conda: '/envs/mapping.yml'
}

process BOWTIE2_INTERLEAVED {
    //
}

process SAMTOOLS_SORT {
    //
}

process SAMTOOLS_INDEX {
    //
}

process CALL_VARIANTS {
    //
}

process FILTER_VARIANTS {
    //
}

process INDEX_VARIANTS {
    //
}

process GET_CONSENSUS {
    //
}

process EXTRACT_CONSENSUS {
    //
}
