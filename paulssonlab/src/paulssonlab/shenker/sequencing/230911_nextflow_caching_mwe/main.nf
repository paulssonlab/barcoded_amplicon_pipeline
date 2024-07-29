import static functions.*

process PROCESS1 {
    input:
    tuple val(meta), val(something)

    output:
    tuple val(meta), val(something)

    exec:
    println "RUNNING COMPUTATION"
}

process PROCESS2 {
    cache false

    input:
    val(meta)

    output:
    val(meta)

    exec:
    println "YOU WON'T SEE THIS"
    meta["foo"] = "bar"
}

workflow {
    main:
    def samples_in = Channel.of([x: [0]])
    map_call_process2(PROCESS1,
                     samples_in,
                     ["x"],
                     [id: { pod5, meta -> pod5 }],
                     "x",
                     ["y"],
                     "_key") { pod5, meta -> [pod5] }
        | PROCESS2
        | view

}
