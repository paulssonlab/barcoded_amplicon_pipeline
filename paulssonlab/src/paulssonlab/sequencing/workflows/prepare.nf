import java.nio.file.Paths
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

include { scp;
          join_key;
          join_each;
          join_map;
          edit_map_key;
          file_in_dir } from '../functions.nf'

include { ANY2FASTA;
          MERGE_FASTAS } from '../modules/fastas.nf'

workflow PREPARE_SAMPLE_SHEET {
    main:
    def sample_sheet_path = Paths.get(params.data_dir, params.sample_sheet_name)
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
           [it, scp("${remote_path}/${it}", Paths.get(params.data_dir, it))]
        }
        .set { ch_reads }
    join_key(samples_in, ch_reads, "reads_path", "reads")
        .set { samples }

    emit:
    samples
}

def json_command(command, input) {
    def SCRIPT_TIMEOUT = 100000 // msec
    def json_input = JsonOutput.toJson(input)
    def proc = command.execute()
    def output_stream = new StringBuffer();
    proc.consumeProcessOutput(output_stream, System.err)
    proc.withWriter { writer ->
        writer.write(json_input)
    }
    proc.waitForOrKill(SCRIPT_TIMEOUT)
    def json_slurper = new JsonSlurper()
    if (proc.exitValue() != 0) {
        println "${command} failed with output:"
        println output_stream.toString()
        throw new Exception("${command} failed")
    }
    def json_output = json_slurper.parseText(output_stream.toString())
    return json_output
}

def get_registry_seqs(references_dir, sample_references) {
    def output = json_command(["${src}/sequencing/bin/get_registry_seqs.py", "${src}/shenker/cloning", references_dir], sample_references)
    edit_map_key(output, "references", "references") { refs ->
        refs.collect { ref -> file_in_dir(references_dir, ref) } as Set
    }
}

workflow PREPARE_REFERENCES {
    take:
    samples_in

    main:
    def references_dir = Paths.get(params.references_dir)
    samples_in
        .map { [(it.run_path): it.reference_names] }
        .collect()
        .map { it.sum() }
        .map { get_registry_seqs(references_dir, it) }
        .set { ch_get_registry_seqs }

    // UNIQUE REFERENCES (convert to fasta, publishDir)
    ch_get_registry_seqs
        .map { it.values()*.references.sum().unique() }
        .flatMap { it.collect { f -> [[id: f.getBaseName()], f] } }
        .set { ch_references_orig_format }

    ch_references_orig_format
        | ANY2FASTA
        | set { ch_references_fasta }

    ch_references_fasta
        .map { [(it[0].id): it[1]] }
        .collect()
        .map { it.sum() }
        .set { ch_references_fasta_map }

    // UNIQUE REFERENCE_SETS (merge fasta)
    ch_get_registry_seqs
        .flatMap { it.values()*.references.unique()*.collect { f -> f.getBaseName() } }
        .map { [id: it, references: it] }
        .set { ch_reference_sets_orig_format }

    join_each(ch_reference_sets_orig_format, ch_references_fasta_map, "references", "references")
        .map { [[id: it.id], it.references] }
        .set { ch_reference_sets }

    MERGE_FASTAS(ch_reference_sets)
        .set { ch_merged_references }

    ch_merged_references
        .view()

    ch_get_registry_seqs
        .map { refs ->
            edit_map_key(refs, "references", "reference_basenames") { it*.getBaseName() }
        }
        .set { ch_get_registry_seqs_with_basenames }

    join_map(samples_in, ch_get_registry_seqs_with_basenames, "run_path")
        .set { samples }

    emit:
    samples
    references = ch_merged_references
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
