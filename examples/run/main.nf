include { MAIN } from "${PIPELINE_ROOT}/sequencing/workflows/nanopore_fish/main.nf"

workflow {
    MAIN()
}
