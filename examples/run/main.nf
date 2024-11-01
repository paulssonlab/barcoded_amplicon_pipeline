include { MAIN } from "${PIPELINE_ROOT}/paulssonlab/src/paulssonlab/sequencing/workflows/nanopore_fish/main.nf"

workflow {
    MAIN()
}
