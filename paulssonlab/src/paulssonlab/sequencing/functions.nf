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

// ch_a should be a channel containing maps of format [old_key_a:value_to_join_on, ...]
// ch_b should be a channel containing lists of format [key, value, ignored...]
// output will be same as ch_a but with old_key renamed to new_key
// and its corresponding value mapped according to the [key, value] lookup table of ch_b
def join_key(ch_a, ch_b, old_key, new_key) {
    ch_b.cross(ch_a.map { [it[old_key], it] } ).map { [*:it[1][1], (new_key): it[0][1]] }
}

// ch_a should be a channel containing maps of format [old_key_a:value_to_join_on, ...]
// same as join_key but ch_b contains lists of format [[old_key_b:value_to_join_on, ...], new_value, ignored...]
// each list is first mapped to a list [k, v] by indexing into id_map with old_key_b
def join_key2(ch_a, ch_b, old_key_a, old_key_b, new_key) {
    join_key(ch_a, ch_b.map { [it[0][old_key_b], it[1]] }, old_key_a, new_key)
}

// each value in ch_a is intepreted as a map of lists of values to join on,
// and each value of that list is mapped according to the [key, value] lookup table of ch_b
def join_each(ch_a, ch_b, old_key, new_key) {
    ch_a.combine(ch_b).map { a, b ->
        [*:a, (new_key): a[old_key].collect { b.get(it) }]
    }
}

// ch_entries is a map of format [key:value_to_join_on, ...]
// ch_map is a channel which contains a single value
// which should be a map [val1:map1, val2:map2, val3:map3, ...]
// the map corresponding to value_to_join_on is merged with the original map from ch_entries
def join_map(ch_entries, ch_map, key) {
    ch_entries.combine(ch_map).map { entry, map ->
        [*:entry, *:map.getOrDefault(entry[key], [:])]
    }
}

// edit_map_key renames a map's key (if closure is Closure.IDENTITY)
// or runs a closure on a map's key (if old_key == new_key)
// or does both at the same time
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

// removes keys from map
def remove_keys(map, keys) {
    def new_map = map.clone()
    new_map.keySet().removeAll(keys)
    new_map
}

// map should be a map [k1:closure1, k2:closure2, ...]
// returns a map [k1:closure1(x), k2:closure2(x), ...]
// for some fixed value x
def collect_closures(map, x) {
    map.collectEntries { k, v -> [(k): v(x) ] }
}

// see below. call_closure is the same as call_process except it does not uniquify the channel
// sent to the process.
def call_closure(Closure closure, ch, join_keys, closure_map, output_keys, Closure preprocess) {
    def ch_input = ch.map { it ->
            def join_map_ = it.subMap(join_keys)
            // we could use preprocess(it) and collect_closures(..., it) here,
            // but that would allow the process to depend on information
            // beyond what's contained in the join keys
            [[*:join_map_,
              *:collect_closures(closure_map, join_map_)], *preprocess(join_map_)]
        }
    def ch_output = closure(ch_input).map { [remove_keys(it[0], closure_map.keySet()), *it[1..-1]] }
    def ch_orig = ch.map { [it.subMap(join_keys), it] }
    ch_output.cross(ch_orig)
        .map {
            [*:it[1][1], *:[output_keys, it[0][1..-1]].transpose().collectEntries()]
        }
}

// preprocess takes each map from ch, removing all keys except for join_keys, and outputs a list
// that's used as input into the process
// a map containing only join_keys will be prepended to this list before sending it to the channel
// process should take this as its first input and pass it through as its first output
// usually it's called "meta"
// the remaining outputs of process will be added to the maps in ch and given the names output_keys
// closure_map is a way of adding additional keys to meta derived from only the information
// in join_keys. the purpose is to allow using processes that expect extra keys in meta
// without having to modify the process itself.
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
