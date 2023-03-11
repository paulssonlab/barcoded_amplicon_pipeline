import java.nio.file.Paths
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

import static functions.*
// include { scp;
//           join_key;
//           join_key2;
//           join_each;
//           join_map;
//           edit_map_key;
//           file_in_dir;
//           map_call_process } from '../functions.nf'

include { ANY2FASTA;
          MERGE_FASTAS } from '../modules/fastas.nf'

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

    join_map(samples_in, ch_get_registry_seqs, "run_path")
        .set { ch_samples_with_references }

    map_call_process(ANY2FASTA,
                     ch_samples_with_references,
                     ["references"],
                     [id: { ref, meta -> ref.baseName }],
                     "references",
                     ["references_fasta"],
                     "_ANY2FASTA") { ref, meta -> [ref] }
        .set { ch_samples_with_fastas }

    call_process(MERGE_FASTAS,
                 ch_samples_with_fastas,
                 "references_fasta",
                 [id: { it.references_fasta }],
                 ["reference"]) { [it.references_fasta] }
        .set { samples }

    emit:
    samples
}
