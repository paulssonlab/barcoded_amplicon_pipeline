import java.nio.file.Paths
import java.nio.file.Files
import java.nio.file.FileSystems
import nextflow.util.CsvParser
import static nextflow.Nextflow.file
import static nextflow.Nextflow.groupKey
import nextflow.Channel
import groovy.transform.Field // SEE: https://stackoverflow.com/a/31301183

@Field static final normal = "\033[0;0m"
@Field static final bold = "\033[0;1m"

static def file_in_dir(dir, filename) {
    file(Paths.get(dir as String, filename as String))
}

static boolean is_dir_empty(dir) {
    try {
        def dir_stream = Files.newDirectoryStream(dir)
        return !dir_stream.iterator().hasNext()
    } catch (Exception e) {
        return true
    }
}

static def download_data(params) {
    def remote_path = Paths.get(params.remote_path_base, params.remote_path)
    def dest = file(params.data_dir)
    if (is_dir_empty(dest)) {
        if (!dest.exists()) {
            dest.mkdirs()
        }
        println "${bold}Downloading data from${normal} ${remote_path}"
        rsync(remote_path, dest)
        println "${bold}Done.${normal}"
    }
    return dest
}

static def rsync(remote_path, dest_path) {
    // SEE: https://stackoverflow.com/questions/25300550/difference-in-collecting-output-of-executing-external-command-in-groovy
    // and https://stackoverflow.com/questions/159148/groovy-executing-shell-commands
    // without escaping spaces, would need -T to deal with scp quoted source paths
    // SEE: https://stackoverflow.com/questions/54598689/scp-fails-with-protocol-error-filename-does-not-match-request
    // not sure why I need two backslashes to escape spaces
    def command = ['rsync', '-az', (remote_path as String).replaceAll(' ', '\\\\ ') + "/", dest_path]
    def proc = command.execute()
    def outputStream = new StringBuffer()
    proc.waitForProcessOutput(outputStream, outputStream)
    if (proc.exitValue() != 0) {
        println "${bold}rsync failed with output:${normal}"
        println outputStream.toString()
        throw new Exception("rync of '${remote_path}' to '${dest_path}' failed")
    }
}

static def glob(pattern, base) {
    def paths = []
    def path_matcher = FileSystems.getDefault().getPathMatcher("glob:" + pattern)
    base = file(base).toAbsolutePath()
    Files.walk(base).forEach { abs_path ->
        def rel_path = base.relativize(abs_path)
        if (path_matcher.matches(rel_path)) {
            paths << abs_path
        }
    }
    return paths.sort() // useful for making nextflow hashes insensitive to input order
}

static def glob_inputs(ch, base, input_names) {
    ch.map { it ->
        input_names.each { k ->
            it[k] = (it.get(k) ? glob(it[k], base) : [])
        }
        it
    }
}

static def find_inputs(ch, base, input_names) {
    ch.map { it ->
        input_names.each { k ->
            if (it.get(k) && file_in_dir(base, it[k]).exists()) {
                it[k] = file_in_dir(base, k)
            } else {
                it.remove(k)
            }
        }
        it
    }
}

static def chunk_files(files, chunk_bytes = 0, chunk_files = 0) {
    if (chunk_bytes && !chunk_files) {
        def chunks = [[]]
        def bytes = 0
        files.each {
            def file_bytes = it.size()
            if (bytes + file_bytes >= chunk_bytes) {
                chunks << []
                bytes = 0
            }
            bytes += file_bytes
            chunks[-1] << it
        }
        chunks
    } else if (chunk_files && !chunk_bytes) {
        files.collate(chunk_files as int)
    } else {
        throw new Exception("exactly one of chunk_bytes, chunk_files must be specified")
    }
}

static def read_tsv(path) {
    def parser = new CsvParser()
        .setSeparator('\t')
    def data = []
    def tsv_file = file(path)
    tsv_file.eachLine { line ->
        data << parser.parse(line)
    }
    return data
}

static def get_samples(params, defaults = [:], substitute = true, Closure preprocess = null) {
    def sample_list = SampleSheetParser.load(params.samples, defaults, substitute, preprocess)
    return Channel.fromList(sample_list).map {
        [*:it, output_run_dir: file_in_dir(params.output, it.run_path)]
    }
}

// equivalent Groovy's GStringImpl and Java's String do not hash to the same value
// so mismatches can cause problems be used as keys in HashMaps, cross, join_key, etc.
// this function converts all GStringImpl's to String and leaves all other objects unchanged
static def stringify(str) {
    (str instanceof GString) ? str.toString() : str
}

// equivalent to ch_a.cross(ch_b) but optionally allows stringifying keys
static def cross(ch_a, ch_b, stringify_keys = true) {
    if (stringify_keys) {
        ch_a = ch_a.map { [stringify(it[0]), *it[1..-1]] }
        ch_b = ch_b.map { [stringify(it[0]), *it[1..-1]] }
    }
    ch_a.cross(ch_b)
}

// ch_a should be a channel containing maps of format [old_key_a:value_to_join_on, ...]
// ch_b should be a channel containing lists of format [key, value, ignored...]
// output will be same as ch_a but with old_key renamed to new_key
// and its corresponding value mapped according to the [key, value] lookup table of ch_b
static def join_key(ch_a, ch_b, old_key, new_key, stringify_keys = true) {
    cross(ch_b, ch_a.map { [it[old_key], it] }, stringify_keys).map { [*:it[1][1], (new_key): it[0][1]] }
}

// ch_a should be a channel containing maps of format [old_key_a:value_to_join_on, ...]
// same as join_key but ch_b contains lists of format [[old_key_b:value_to_join_on, ...], new_value, ignored...]
// each list is first mapped to a list [k, v] by indexing into id_map with old_key_b
static def join_key2(ch_a, ch_b, old_key_a, old_key_b, new_key, stringify_keys = true) {
    join_key(ch_a, ch_b.map { [it[0][old_key_b], it[1]] }, old_key_a, new_key, stringify_keys)
}

// each value in ch_a is intepreted as a map of lists of values to join on,
// and each value of that list is mapped according to the [key, value] lookup table of ch_b
static def join_each(ch_a, ch_b, old_key, new_key, stringify_keys = true) {
    if (stringify_keys) {
        ch_b = ch_b.map { it.collectEntries { k, v -> [(stringify(k)): v] } }
    }
    ch_a.combine(ch_b).map { a, b ->
        [*:a, (new_key): a[old_key].collect { b.get(stringify_keys ? stringify(it) : it) }]
    }
}

// ch_entries is a map of format [key:value_to_join_on, ...]
// ch_map is a channel which contains a single value
// which should be a map [val1:map1, val2:map2, val3:map3, ...]
// the map corresponding to value_to_join_on is merged with the original map from ch_entries
static def join_map(ch_entries, ch_map, key, stringify_keys = true) {
    if (stringify_keys) {
        ch_map = ch_map.map { it.collectEntries { k, v -> [(stringify(k)): v] } }
    }
    ch_entries.combine(ch_map).map { entry, map ->
        [*:entry, *:map.getOrDefault(stringify_keys ? stringify(entry[key]) : entry[key], [:])]
    }
}

// returns a map with old_key removed and replaced by new_key, optionally
// modifying the value with a closure
static def rename_key(map, old_key, new_key, Closure closure = Closure.IDENTITY) {
    def new_map = map.clone()
    new_map.keySet().removeAll([old_key])
    new_map[new_key] = closure(map[old_key])
    new_map
}

static def swap_key(map, key1, key2, key3 = null) {
    key3 = key3 ?: key1
    def temp = map[key1]
    def new_map = map.close()
    new_map.put(key1, new_map[key2]);
    new_map.put(key2, temp);
    new_map
}

static def edit_key(map, key, Closure closure = Closure.IDENTITY) {
    [*:map, (key): closure(map.get(key))]
}

// removes keys from map
static def remove_keys(map, keys) {
    def new_map = map.clone()
    new_map.keySet().removeAll(keys)
    new_map
}

// map should be a map [k1:closure1, k2:closure2, ...]
// returns a map [k1:closure1(x), k2:closure2(x), ...]
// for some fixed value x
static def collect_closures(map, Object... args) {
    map.collectEntries { k, v -> [(k): v(*args)] }
}

// see below. call_closure is the same as call_process except it does not uniquify the channel
// sent to the process.
static def call_closure(Closure closure, ch, join_keys, closure_map, output_keys, when = null, passthrough = true, stringify_keys = true, Closure preprocess) {
// static def call_closure(Closure closure, ch, join_keys, closure_map, output_keys, Closure when, passthrough = true, stringify_keys = true, Closure preprocess) {
    if (!(output_keys instanceof List)) {
        throw new Exception("output_keys must be a List")
    }
    def ch_filtered
    def ch_passthrough
    if (when) {
        ch_filtered = ch.filter { when(it) }
        ch_passthrough = ch.filter { !when(it) }
    } else {
        ch_filtered = ch
        ch_passthrough = Channel.empty()
    }
    def ch_input = ch_filtered.map { it ->
            def join_map_ = it.subMap(join_keys)
            // we could use preprocess(it) and collect_closures(..., it) here,
            // but that would allow the process to depend on information
            // beyond what's contained in the join keys
            [[*:join_map_,
              *:collect_closures(closure_map, join_map_)], *preprocess(join_map_)]
        }
    def ch_processed = closure(ch_input).map { [remove_keys(it[0], closure_map.keySet()), *it[1..-1]] }
    def ch_orig = ch.map { [it.subMap(join_keys), it] }
    return cross(ch_processed, ch_orig, stringify_keys)
        .map {
            [*:it[1][1], *:[output_keys, it[0][1..-1]].transpose().collectEntries()]
        }
        .mix(ch_passthrough)
}

// preprocess takes each map from ch (with all keys removed except for join_keys) and outputs a list
// that's used as input into the process
// a map containing only join_keys will be prepended to this list before sending it to the channel
// process should take this as its first input and pass it through as its first output
// usually it's called "meta"
// the remaining outputs of process will be added to the maps in ch and given the names output_keys
// closure_map is a way of adding additional keys to meta derived from only the information
// in join_keys. the purpose is to allow using processes that expect extra keys in meta
// without having to modify the process itself.
static def call_process(process, ch, join_keys, closure_map, output_keys, when = null, passthrough = true, stringify_keys = true, Closure preprocess) {
// static def call_process(process, ch, join_keys, closure_map, output_keys, Closure when, passthrough = true, stringify_keys = true, Closure preprocess) {
    // TODO: not sure why { process(it.unique()) } doesn't work
    // call_closure({ it.unique() | process }, ch, join_keys, closure_map, output_keys, when, passthrough, stringify_keys, preprocess)
    call_closure({ it.unique() | process }, ch, join_keys, closure_map, output_keys, when, passthrough, stringify_keys, preprocess)
}

// this works like call_process but for each input map emitted by ch, it will call process
// for each value in a collection (i.e., list) in the input map under map_input_key
// all arguments work similarly to those for call_process, the only difference is that
// closure_map and preprocess are allowed to depend on the mapped collection
// and not just the mapped value itself, although this may prevent deduplicating process calls
// temp_key is an arbitrary unique string and is used as a prefix for keys that are
// temporarily added to the map that is passed through process
static def map_call_process(process, ch, join_keys, closure_map, map_input_key, output_keys, temp_key, stringify_keys = true, passthrough = true, Closure preprocess) {
    if (!(output_keys instanceof List)) {
        throw new Exception("output_keys must be a List")
    }
    // the key we use to store the list of hashes we use when joining input channel and process output channel
    def temp_hashes_key = "${temp_key}_hashes"
    // the key we use to store mapped value for each process input. this is used to ensure
    // that the order of values in output collections corresponds to the order of values
    // in the input collection
    def temp_mapped_value_key = "${temp_key}_mapped_value"
    // ch is a channel emitting maps. some of these items may contains empty collections.
    // these items disappear when we transpose, below. as such, we filter them out.
    // if passthrough is true, we mix them to ch_output before returning
    def ch_input_nonempty = ch.filter { it.getOrDefault(map_input_key, []).size() != 0 }
    def ch_input_untransposed = ch_input_nonempty.map {
            def collection_to_map = it.getOrDefault(map_input_key, [])
            // we need to use groupKey to wrap the hash so that the second groupTuple invokation
            // (ch_output_transposed.groupTuple) knows how many elements to expect for each hash
            // note that Nextflow's transpose operator requires that the collection be a List
            [groupKey(it.hashCode(), collection_to_map.size()), it, collection_to_map as List]
        }
    // transpose over collection (all tuples for a given collection get the same hash)
    def ch_input_transposed = ch_input_untransposed.transpose(by: 2)
    def ch_input_transposed_with_key = ch_input_transposed.map { key, map, mapped_value ->
        // this is the map that is passed as the first argument to the process
        // (additional keys are added according to closure_map)
        def join_map_ = map.subMap(join_keys)
        // this is the map that is passed to the closures (preprocess and values of closure_map)
        def closure_input_map = map.subMap([*join_keys, map_input_key])
        // note that the closures (preprocess and values of closure_map) can in principle
        // depend on the collection we're mapping over (specified by map_input_key)
        // but that these dependencies may prevent deduplicating process calls
        // process_map, which is the first item in the tuple passed as input to the process,
        // includes the values in join_map_, the output of collect_closures,
        // as well as mapped_value
        def process_map = [*:join_map_, *:collect_closures(closure_map, mapped_value, closure_input_map),
                           (temp_mapped_value_key):mapped_value]
        // preprocess is a closure that returns a tuple of all additional inputs to the process
        // besides process_map
        def process_input = [process_map, *preprocess(mapped_value, closure_input_map)]
        [process_input, key]
    }
    // we want to get a list of hashes for each unique process_input
    // (so that we can reuse output for identical computations)
    def ch_process_groups = ch_input_transposed_with_key.groupTuple(by: 0)
    def ch_process_input = ch_process_groups.map { process_input, keys ->
        // we have to store the list of hashes in the map so that we have it in the process output
        [[*:process_input[0], (temp_hashes_key):keys], *process_input[1..-1]]
    }
    def ch_process_output = ch_process_input | process
    // we need to pull out the list of hashes and put it in a tuple so we can transpose
    def ch_output_untransposed = ch_process_output.map {
        [it[0].get(temp_hashes_key), it]
    }
    def ch_output_transposed = ch_output_untransposed.transpose(by: 0)
    // groupTuple will emit each output collection
    // because we are using groupTuple with groupKeys that have their sizes specified, this channel
    // will emit as soon as all output elements for an input collection are ready
    def ch_output_grouped = ch_output_transposed.groupTuple(by: 0)
    // join the inputs (with hashes) and the outputs
    def ch_output_joined = ch_input_untransposed.cross(ch_output_grouped).map { input, output ->
        // remove hashes from input/output tuples
        [input[1], output[1]]
    }
    def ch_output = ch_output_joined.map { input, output ->
        // make a map from input_value => (output1, output2, ...)
        def outputs_by_mapped_values = output.collectEntries { [(it[0].get(temp_mapped_value_key)): it[1..-1]] }
        // make output collections ordered according to the input collection
        def new_output = output_keys.withIndex().collectEntries { output_key, idx ->
            [(output_key): input.get(map_input_key).collect { mapped_value ->
                outputs_by_mapped_values.get(mapped_value)[idx]
                }]
        }
        // merge new output keys into map
        [*:input, *:new_output]
    }
    if (passthrough) {
        // add items with empty collections to output channel
        // we could add an empty list to each item keyed under output_key,
        // but this probably is not helpful
        def ch_input_empty = ch.filter { it.getOrDefault(map_input_key, []).size() == 0 }
        ch_output.mix(ch_input_empty)
    } else {
        ch_output
    }
}
