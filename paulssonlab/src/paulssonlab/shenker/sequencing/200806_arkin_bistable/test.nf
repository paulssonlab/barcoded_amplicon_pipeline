nextflow.enable.dsl=2

def join_key(ch_a, ch_b, key) {
    return ch_a
}

def join_each(ch_a, ch_b, old_key, new_key) {
    return ch_a
}

// TODO: what is this useful for??
def inner_join(ch_a, ch_b) {
    return ch_b.cross(ch_a).map { [it[0][0], *it[1][1..-1], *it[0][1..-1]] }
}

workflow A {
    Channel.fromList([[1, "a", "aa"], [2, "b", "bb"], [3, "c", "cc"], [4, "d", "dd"]])
        .set { reads }
    Channel.fromList([[id: "a", reads: 4],
                      [id: "b", reads: 3],
                      [id: "c", reads: 2],
                      [id: "d", reads: 1]])
        .set { samples }
    join_key(samples, reads, "reads").view()
    // expected output: [1, [id: "d", reads: 1], "a", "aa"]
}

workflow B {
    Channel.of([2: "a2", 3: "b3", 4: "c4", 5: "d5", 6: "e6"])
        .set { refs }
    Channel.fromList([[ref_names: [2,3,4]], [ref_names: [4,5,6]]])
        .set { samples }
    join_each(samples, refs, "ref_names", "ref_files").view()
    // expected output: [ref_names: [2,3,4], ref_files: ["a2", "b3", "c4"]]
}

workflow {
    A()
}
