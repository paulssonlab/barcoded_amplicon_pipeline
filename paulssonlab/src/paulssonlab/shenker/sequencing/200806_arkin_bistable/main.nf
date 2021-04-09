nextflow.enable.dsl=2

include { MAIN } from "${src}/sequencing/workflows/illumina_whole_plasmid/main.nf"

workflow {
    MAIN()
}
