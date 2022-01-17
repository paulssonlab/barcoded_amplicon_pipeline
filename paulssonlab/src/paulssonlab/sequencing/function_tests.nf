nextflow.enable.dsl=2

include { join_key; join_key2;
          join_each; join_map;
          edit_map_key; remove_keys;
          collect_closures; call_closure;
          call_process; map_call_process } from './functions.nf'

workflow JOIN_KEY {
    Channel.fromList([[id: "a", reads: 4],
                      [id: "b", reads: 3],
                      [id: "c", reads: 2],
                      [id: "d", reads: 1]])
        .set { samples }
    // only the first two elements of each list are used
    Channel.fromList([[1, "a", "xx"],
                      [2, "b", "xx"],
                      [3, "c", "xx"],
                      [4, "d", "xx"]])
        .set { reads }
    join_key(samples, reads, "reads", "reads_files").view()
    // expected output: [reads_files:a, id:d, reads:1], ...
}

workflow JOIN_KEY2 {
    Channel.fromList([[id: "a", reads: 4],
                      [id: "b", reads: 3],
                      [id: "c", reads: 2],
                      [id: "d", reads: 1]])
        .set { samples }
    Channel.fromList([[[read_id: 1, blah: "x"], "a", "xx"],
                      [[read_id: 2, blah: "x"], "b", "xx"],
                      [[read_id: 3, blah: "x"], "c", "xx"],
                      [[read_id: 4, blah: "x"], "d", "xx"]])
        .set { reads }
    join_key2(samples, reads, "reads", "read_id", "reads_files").view()
    // expected output: [reads_files:a, id:d, reads:1], ...
}

workflow JOIN_EACH {
    Channel.fromList([[ref_names: [2,3,4]], [ref_names: [4,5,6]]])
        .set { samples }
    Channel.of([2: "a2", 3: "b3", 4: "c4", 5: "d5", 6: "e6"])
        .set { refs }
    join_each(samples, refs, "ref_names", "ref_files").view()
    // expected output: [ref_names: [2,3,4], ref_files: ["a2", "b3", "c4"]], ...
}

workflow JOIN_MAP {
    Channel.fromList([[id: "a", flag: 1],
                      [id: "b", flag: 2],
                      [id: "c", flag: 3],
                      [id: "d", flag: 4]])
        .set { samples }
    Channel.of([1: [config: "x"], 2: [config: "y"], 3: [config: "z"], 4: [config: "zz"]])
        .set { configs }
    join_map(samples, configs, "flag").view()
    // expected output: [id:a, flag:1, config:x], ...
}

process DUMMY_PROCESS {
    tag "$meta.id"

    input:
    tuple val(meta), val(input1), val(input2)

    output:
    tuple val(meta), path('x'), path('y')

    script:
    """
    echo x input1:${input1} + input2:${input2} + args:${meta.args} + id:${meta.id} > x
    echo y input1:${input1} + input2:${input2} + args:${meta.args} + id:${meta.id} > y
    """
}

workflow CALL_PROCESS {
    ch = Channel.fromList([[val: 1, input1: "foo", input2: "a", args: "foo2", blah: 5],
                           [val: 2, input1: "foo", input2: "a", args: "foo2", blah: 6],
                           [val: 3, input1: "bar", input2: "b", args: "bar2", blah: 7],
                           [val: 4, input1: "bar", input2: "c", args: "bar2", blah: 8]])
    // args: process, ch, join_keys, closure_map, output_keys, Closure preprocess
    call_process(DUMMY_PROCESS,
                ch,
                ["input1", "input2", "args"],
                [id: { it.val }],
                ["out1", "out2"]) { [it.input1, it.input2] }
        .view()
    // expected output:
    // [val:1, input1:foo, input2:a, args:foo2, blah:5, out1:/tmp/nextflow/ad/3ac3e469f4889b15a3c71af8a8f390/x, out2:/tmp/nextflow/ad/3ac3e469f4889b15a3c71af8a8f390/y]
    // [val:2, input1:foo, input2:a, args:foo2, blah:6, out1:/tmp/nextflow/ad/3ac3e469f4889b15a3c71af8a8f390/x, out2:/tmp/nextflow/ad/3ac3e469f4889b15a3c71af8a8f390/y]
    // [val:3, input1:bar, input2:b, args:bar2, blah:7, out1:/tmp/nextflow/87/64875d579ea44aa6506b29b6b3f0fc/x, out2:/tmp/nextflow/87/64875d579ea44aa6506b29b6b3f0fc/y]
    // [val:4, input1:bar, input2:c, args:bar2, blah:8, out1:/tmp/nextflow/df/d874813693851e31d3d42f117e374b/x, out2:/tmp/nextflow/df/d874813693851e31d3d42f117e374b/y]
}

//def map_call_process(process, ch, join_keys, closure_map, map_input_keys, map_output_keys, Closure preprocess) {
//    return 0
    // def closure = {
    //     it.map {println x; x}
    //     // it
    //     // it.unique() | process
    //     // it.map {
    //     //     []
    //     // }
    // }
    // call_closure(closure, ch, join_keys, closure_map, map_output_keys, Closure.IDENTITY)
//}

process DUMMY_PROCESS_INDEX {
    tag "$meta.id"

    input:
    tuple val(meta), val(input1), val(input2)

    output:
    tuple val(meta), path('x'), path('y')

    script:
    """
    echo x input1:${input1} + input2:${input2} + args:${meta.args} + id:${meta.id} > x
    echo y input1:${input1} + input2:${input2} + args:${meta.args} + id:${meta.id} > y
    """
}

workflow MAP_CALL_PROCESS {
    ch = Channel.fromList([[val: 1, args: "foo", refs1: ["k", "kk", "kkk"], refs2: ["a", "aa", "aaa"]],
                           [val: 2, args: "foo", refs1: ["k", "kk", "kkk"], refs2: ["a", "aa", "aaa"]],
                           [val: 3, args: "bar", refs1: ["k", "kk", "kkk"], refs2: ["a", "aa", "aaa"]],
                           [val: 4, args: "bar", refs1: ["l", "ll", "lll"], refs2: ["a", "aa", "aaa"]]])
    // args: process, ch, join_keys, closure_map, map_input_key, map_output_keys, Closure preprocess
    map_call_process(DUMMY_PROCESS,
                     ch,
                     ["args"],
                     [id: { it.refs }],
                     ["refs1", "refs2"],
                     ["out1", "out2"]) { meta -> [meta.refs1, meta.refs2] }
        .view()
}
