nextflow.enable.dsl=2

@Grab(group='com.moandjiezana.toml', module='toml4j', version='0.7.2')
import com.moandjiezana.toml.Toml;

include { MAIN } from "${src}/sequencing/workflows/illumina_whole_plasmid/main.nf"

process MAKE_INDEX {
    input:
    tuple val(meta)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    shell:
    '''
    echo !{meta.id} > !{bam}
    '''
}

workflow {
    //MAKE_INDEX()
}
