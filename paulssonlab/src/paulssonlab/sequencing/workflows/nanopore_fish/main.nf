import java.nio.file.Files
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import static functions.*

workflow NANOPORE_FISH {
    take:
    samples_in

    main:
    samples_in
        | map {
            if (it.get("filtered_reads")) {
                rename_key(rename_key(it, "reads", "unfiltered_reads"), "filtered_reads", "reads")
            } else {
                it
            }
        }
        | call_MINIMAP2_ALIGN
        | call_SAMTOOLS_INDEX
        | map {
            it.output_run_dir.mkdirs()
            ["reference", "reads", "bam", "bam_index"].each { k ->
                (it[k] instanceof Collection ? it[k] : [it[k]]).each { src ->
                    Files.copy(src, file_in_dir(it.output_run_dir, src.name),
                               StandardCopyOption.REPLACE_EXISTING)
                }
            }
            it
        }
        | set { samples }
    // pod5 merging (if necessary)
    // pod5 splitting
    // dorado simplex+duplex model download
    // dorado duplex basecalling
    // publish bams
    // filter dx:0/1 from bam -> fastq
    // porechop_abi
    // minigraph barcode_v1.gfa

    samples.view()

    emit:
    samples

}

workflow MAIN {
    // download_data(params) // TODO: remove?
    samples_preglob = get_samples(params, [:], true)
    s2 = glob_inputs(samples_preglob, params.root, ["fastq_input", "fast5_input", "pod5_input"])
    find_inputs(s2, params.root, ["fastq_input2", "fast5_input2", "pod5_input2"])
        // | NANOPORE_FISH
        | view
        | set { samples }

    emit:
    samples
}

workflow {
    MAIN()
}
