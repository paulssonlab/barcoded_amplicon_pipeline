def join_each(ch_a, ch_b):
    return 0

def join_key(ch_a, ch_b, key):
    return 0

// TODO: what is this useful for??
def inner_join(ch_a, ch_b) {
    return ch_b.cross(ch_a).map { [it[0][0], *it[1][1..-1], *it[0][1..-1]] }
}

workflow A {
    Channel.of([2: "a2", 3: "b3", 4: "c4", 5: "d5", 6: "e6"])
        .set { refs }
    Channel.fromList([[refs: [2,3,4]], [refs: [4,5,6]]]])
        .set { samples }
    join_each(samples, refs, "refs").view()
}

workflow B {
    Channel.fromList([[1, "a"], [2, "b"], [3, "c"], [4, "d"]])
        .set { reads }
    Channel.fromList([[id: "a", reads: 4],
                      [id: "b", reads: 3],
                      [id: "c", reads: 2],
                      [id: "d", reads: 1]])
        .set { samples }
    join_key(samples, reads, "reads")
}

workflow {
    B()
}
