include { read_tsv } from '../functions.nf'

include { ANY2FASTA;
          MERGE_FASTAS } from '../modules/fastas.nf'

workflow PREPARE_SAMPLE_SHEET {
    main:
    def sample_sheet_path = "${params.data_dir}/${params.sample_sheet_name}"
    def remote_path = "${params.remote_path_base}/${params.remote_path}"
    scp("${remote_path}/${params.sample_sheet_name}", sample_sheet_path)
    def samples = read_tsv(sample_sheet_path)

    emit:
    samples
}

workflow PREPARE_READS {
    main:
    // TODO: concatenate reads lists if reads is a list (using flatten)
    // TODO: uniquify reads
    // TODO: scp
    def remote_path = "${params.remote_path_base}/${params.remote_path}"
    for (row in samples) {
        def filename = "${row[0]}.fastq"
        scp("${remote_path}/${filename}", "${params.data_dir}/${filename}")
    }
    // TODO: merge turn reads [DIR OR FILE! possibly a list!] into file() objects
    //samples = Channel.value(sample_list)
    // sample_list
    //     //.splitCsv(sep:'\t')
    //     .map { row -> [id: row[0], references: row[1].split('\s*,\s*') as Set] }
    //     .map { [reads: file("${data_dir}/${it.id}.fastq"), *:it]}
    //     .set { ch_samples }
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
