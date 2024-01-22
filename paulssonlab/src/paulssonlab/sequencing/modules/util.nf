process CAT {
    tag "$meta.id"
    //label "local"

    time 10.min
    memory 200.MB

    input:
    tuple val(meta), path(input, stageAs: "input/*")

    output:
    tuple val(meta), path("${meta.id}")

    script:
    """
    cat ${input} > ${meta.id}
    """

    stub:
    """
    cat ${input} > ${meta.id}
    """
}
