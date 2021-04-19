import java.nio.file.Paths
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

include { scp;
          join_key;
          join_each;
          join_map;
          collect_map_key } from '../functions.nf'

include { ANY2FASTA;
          MERGE_FASTAS } from '../modules/fastas.nf'

workflow PREPARE_SAMPLE_SHEET {
    main:
    def sample_sheet_path = Paths.get(workDir as String, params.data_dir, params.sample_sheet_name)
    def remote_path = Paths.get(params.remote_path_base, params.remote_path)
    scp("${remote_path}/${params.sample_sheet_name}", sample_sheet_path)
    def sample_list = SampleSheetParser.load(sample_sheet_path as String)
    Channel
        .fromList(sample_list)
        .set { samples }

    emit:
    samples
}

workflow PREPARE_READS {
    take:
    samples_in

    main:
    def remote_path = Paths.get(params.remote_path_base, params.remote_path)
    samples_in
        .map { it.reads_path }
        .unique()
        .map {
           [it, scp("${remote_path}/${it}", Paths.get(workDir as String, params.data_dir, it))]
        }
        .set { ch_reads }
    join_key(samples_in, ch_reads, "reads_path", "reads")
        .set { samples }

    emit:
    samples
}

def json_command(command, input) {
    def json_input = JsonOutput.toJson(input)
    def jsonSlurper = new JsonSlurper()
    def str = '''{"default/4309_APA4309_321712w_AE5":{"reference_names":["pLIB219","pLIB220"],"references": ["aaa.gb", "bbb.gb"]},
                  "default/4310_APA4310_321712w_AE6":{"reference_names":["pLIB219","pLIB222"],"references": ["xaa.gb", "bbb.gb"]},
                  "default/4311_APA4311_321712w_AE7":{"reference_names":["pLIB221","pLIB222"],"references": ["yaa.gb", "bbb.gb"]}}'''
    def json_output = jsonSlurper.parseText(str)
    return json_output
}

workflow PREPARE_REFERENCES {
    take:
    samples_in

    main:
    def data_dir = Paths.get(workDir as String, params.data_dir)
    // samples_in.view()
    samples_in
        // .flatMap { [(it.run_path): it.reference_names] }
        .map { [(it.run_path): it.reference_names] }
        .collect()
        .map { it.sum() }
        .map {
            def refs = json_command(["${src}/bin/get_registry_seqs.py", data_dir], it)
            collect_map_key(refs, "references") { ref ->
                file(Paths.get(data_dir as String, ref))
            }
        }
        // .view()
        .set { ch_registry_seqs }
    join_map(samples_in, ch_registry_seqs, "run_path")
        .view()

    // ch_samples
    //     .flatMap { meta -> meta.references }
    //     .unique()
    //     .collect()
    //     .set { ch_reference_names }

    // GET_REGISTRY_SEQS(ch_reference_names.map { it.join(",") })
    //     .flatMap { it.collect { f -> [[id: f.getBaseName()], f] }}
    //     .set { ch_references_orig_format }

    // ANY2FASTA(ch_references_orig_format)
    //     .set { ch_references }

    // ch_references
    //     .map { [[it[0].id, it[1]]] }
    //     .collect()
    //     .map { it.collectEntries() }
    //     .set { ch_collected_references }

    // ch_samples
    //     .map { it.references }
    //     .unique()
    //     .set { ch_reference_sets }

    // ch_reference_sets
    //     .combine(ch_collected_references)
    //     .map { ref_set, refs -> [ref_set, ref_set.collect { refs[it] }] }
    //     .set { ch_reference_set_fastas }

    // MERGE_FASTAS(ch_reference_set_fastas)
    //     .set { ch_merged_references }

    // emit:
    // samples = ch_samples
}

workflow MERGE_INDICES {
    take:
    ch_samples
    ch_indexes

    main:
    inner_join(ch_samples.map { [it.references, it] }, ch_indexes)
        .join(ch_indexes, remainder: true)
        .map { [index: it[2], *:it[1]] }
        .set { ch_samples_indexed }

    emit:
    samples = ch_samples_indexed
    indexes = ch_indexes
}
