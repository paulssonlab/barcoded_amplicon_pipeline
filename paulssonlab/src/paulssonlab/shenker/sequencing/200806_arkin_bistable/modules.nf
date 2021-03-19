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
