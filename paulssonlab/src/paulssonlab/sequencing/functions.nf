import java.nio.file.Paths
import nextflow.util.CsvParser

def file_in_dir(dir, filename) {
    file(Paths.get(dir as String, filename as String))
}

def scp(remote_path, dest_path) {
    def dest = file(dest_path)
    if (!dest.exists()) {
        def dest_parent = dest.getParent()
        if (!dest_parent.exists()) {
            dest_parent.mkdirs()
        }
        // SEE: https://stackoverflow.com/questions/25300550/difference-in-collecting-output-of-executing-external-command-in-groovy
        // and https://stackoverflow.com/questions/159148/groovy-executing-shell-commands
        // without escaping spaces, would need -T to deal with scp quoted source paths
        // SEE: https://stackoverflow.com/questions/54598689/scp-fails-with-protocol-error-filename-does-not-match-request
        // not sure why I need two backslashes to escape spaces
        def command = ['scp', remote_path.replaceAll(' ', '\\\\ '), dest_path]
        def proc = command.execute()
        //def proc = ['scp', "\"${remote_path}\"", dest_path].execute()
        def outputStream = new StringBuffer();
        proc.waitForProcessOutput(outputStream, System.err)
        //proc.waitForProcessOutput(System.out, System.err)
        //return outputStream.toString()
        if ((proc.exitValue() != 0) || (!dest.exists())) {
            println "scp failed with output:"
            println outputStream.toString()
            throw new Exception("scp of '${remote_path}' to '${dest_path}' failed")
        }
    }
    return dest
}

def read_tsv(path) {
    def parser = new CsvParser()
        .setSeparator('\t')
    def data = []
    def tsv_file = file(path)
    tsv_file.eachLine { line ->
        data << parser.parse(line)
    }
    return data
}

def join_key(ch_a, ch_b, old_key, new_key) {
    ch_b.cross(ch_a.map { [it[old_key], it] } ).map { [*:it[1][1], (new_key): it[0][1]] }
}

def join_key2(ch_a, ch_b, old_key_a, old_key_b, new_key) {
    join_key(ch_a, ch_b.map { [it[0][old_key_b], it[1]] }, old_key_a, new_key)
}

def join_each(ch_a, ch_b, old_key, new_key) {
    ch_a.combine(ch_b).map { a, b ->
        [*:a, (new_key): a[old_key].collect { b.get(it) }]
    }
}

def join_map(ch_entries, ch_map, key) {
    ch_entries.combine(ch_map).map { entry, map ->
        [*:entry, *:map.getOrDefault(entry[key], [:])]
    }
}

// TODO: if I define a default argument Closure closure = Closure.IDENTITY,
// I get a wrong number of arguments error
// this may be a limitation of nextflow functions (is there a difference?),
// and real Groovy functions (defined in a .groovy file)
// and/or typed functions may work
def edit_map_key(map, old_key, new_key, Closure closure) {
    map.collectEntries { k, v ->
        def value = [*:v, (new_key): closure(v.get(old_key))]
        [(k): value]
    }
}

def remove_keys(map, keys) {
    def new_map = map.clone()
    new_map.keySet().removeAll(keys)
    new_map
}

def collect_closures(map, x) {
    map.collectEntries { k, v -> [(k): v(x) ] }
}

def call_closure(Closure closure, ch, join_keys, closure_map, output_keys, Closure preprocess) {
    def ch_input = ch.map { it ->
            def join_map = it.subMap(join_keys)
            // we could use preprocess(it) and collect_closures(..., it) here,
            // but that would allow the process to depend on information
            // beyond what's contained in the join keys
            [[*:join_map,
              *:collect_closures(closure_map, join_map)], *preprocess(join_map)]
        }
        //.unique()
        //| process
    def ch_output = closure(ch_input).map { [remove_keys(it[0], closure_map.keySet()), *it[1..-1]] }
    def ch_orig = ch.map { [it.subMap(join_keys), it] }
    ch_output.cross(ch_orig)
        .map {
            [*:it[1][1], *:[output_keys, it[0][1..-1]].transpose().collectEntries()]
        }
}

def call_process(process, ch, join_keys, closure_map, output_keys, Closure preprocess) {
    // TODO: not sure why { process(it.unique()) } doesn't work
    call_closure({ it.unique() | process }, ch, join_keys, closure_map, output_keys, preprocess)
}

// map_call_closure?

def map_call_process(process, ch, join_keys, closure_map, map_input_key, map_output_keys, Closure preprocess) {
    def closure = {
        it.map {println x; x}
        // it
        // it.unique() | process
        // it.map {
        //     []
        // }
    }
    call_closure(closure, ch, join_keys, closure_map, map_output_keys, Closure.IDENTITY)
}

// [meta]
// [[meta], [...,ref2,...]]
// .transpose()
// [[meta], ref1]
// preprocess()
// [[meta], ref1, bowtie2_build_args]
