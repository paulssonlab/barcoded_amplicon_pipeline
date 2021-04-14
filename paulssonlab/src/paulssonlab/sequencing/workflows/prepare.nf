import java.nio.file.Paths

include { scp;
          join_key;
          join_each } from '../functions.nf'

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
    samples

    main:
    def remote_path = Paths.get(params.remote_path_base, params.remote_path)
    samples
        .map { it.reads_path }
        .unique()
        .map {
           [it, scp("${remote_path}/${it}", Paths.get(workDir as String, params.data_dir, it))]
        }
        .set { ch_reads }
    join_key(samples, ch_reads, "reads_path", "reads")
        .set { samples_with_reads }
    // inner_join(samples.map { [it.reads_path, it] }, ch_reads)
    //     .map {}

    // COLLECT, make map
    // MAP OVER reads

    // MAKE FUNCTIONS FOR one-to-one channel joining
    // AND many-to-one

        // .set { read_files }
    // inner_join(samples.map { [it.reads_path, it] }, read_files)
    //     .map { [reads:*:it[1]] }
    // for (row in samples) {
    //     def filename = "${row[0]}.fastq"
    //     scp("${remote_path}/${filename}", "${params.data_dir}/${filename}")
    // }
    // TODO: merge turn reads [DIR OR FILE! possibly a list!] into file() objects
    //samples = Channel.value(sample_list)
    // sample_list
    //     //.splitCsv(sep:'\t')
    //     .map { row -> [id: row[0], references: row[1].split('\s*,\s*') as Set] }
    //     .map { [reads: file("${data_dir}/${it.id}.fastq"), *:it]}
    //     .set { ch_samples }
    emit:
    samples_with_reads
}

workflow PREPARE_REFERENCES {
    take:
    sample_list

    main:

    ch_samples
        .flatMap { meta -> meta.references }
        .unique()
        .collect()
        .set { ch_reference_names }

    GET_REGISTRY_SEQS(ch_reference_names.map { it.join(",") })
        .flatMap { it.collect { f -> [[id: f.getBaseName()], f] }}
        .set { ch_references_orig_format }

    ANY2FASTA(ch_references_orig_format)
        .set { ch_references }

    ch_references
        .map { [[it[0].id, it[1]]] }
        .collect()
        .map { it.collectEntries() }
        .set { ch_collected_references }

    ch_samples
        .map { it.references }
        .unique()
        .set { ch_reference_sets }

    ch_reference_sets
        .combine(ch_collected_references)
        .map { ref_set, refs -> [ref_set, ref_set.collect { refs[it] }] }
        .set { ch_reference_set_fastas }

    MERGE_FASTAS(ch_reference_set_fastas)
        .set { ch_merged_references }

    emit:
    samples = ch_samples
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
