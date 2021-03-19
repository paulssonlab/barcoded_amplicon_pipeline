process SCP {
    input:
    val(remote_path)

    output:
    path('samples.tsv')

    shell:
    '''
    scp !{remote_path} !{}
    '''
}

process GET_REGISTRY_SEQS {
    input:
    stdin

    output:
    path('*.gb', type: 'file')

    storeDir 'test/references'

    shell:
    '''
    #!/usr/bin/env python
    import sys
    import shutil
    import re

    SOURCE = "/Users/jacob/Dropbox (Personal)/Research/Paulsson/paulssonlab/paulssonlab/src/paulssonlab/shenker/sequencing/200806_arkin_bistable/references/pLIB219.gb"

    ids = re.split(r"\\s*,\\s*", sys.stdin.read().rstrip())
    for i, id in enumerate(ids):
        shutil.copy(SOURCE, f"{id}.gb")
        #with open(f"{id}.gb", "w") as f:
        #    f.write(f"BLAH {id} BLEE\\n")
    '''
}

process ANY2FASTA {
    tag "$meta.id"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${meta.id}.fasta")

    publishDir 'test/references'

    conda 'envs/any2fasta.yml'

    shell:
    """
    any2fasta -q $input | seqkit replace -p '(.*)' -r '${meta.id}' > ${meta.id}.fasta
    """
}

process MERGE_FASTAS {
    tag "$meta"

    input:
    tuple val(meta), path('seq')

    output:
    tuple val(meta), path('reference.fasta')

    shell:
    """
    cat seq* > reference.fasta
    """
}

process BOWTIE2_BUILD {
    tag "$meta"

    conda 'envs/mapping.yml'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('bowtie2'), emit: index

    script:
    """
    mkdir bowtie2
    bowtie2-build --threads $task.cpus $fasta bowtie2/${fasta.baseName}
    """
}

process BOWTIE2_INTERLEAVED {
    tag "$meta.id"

    conda 'envs/mapping.yml'

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path('*.bam'), emit: bam
    tuple val(meta), path('*.log'), emit: log

    shell:
    '''
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
    (bowtie2 \
        --threads !{task.cpus} \
        -x $INDEX --interleaved !{reads} \
        | samtools view -@ !{task.cpus} -Sbh -o !{meta.id}.bam -) 2> !{meta.id}.bowtie.log
    '''
}

process SAMTOOLS_SORT {
    tag "$meta.id"

    conda 'envs/mapping.yml'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam

    script:
    """
    samtools sort -@ $task.cpus -o ${meta.id}.bam -T ${meta.id} $bam
    """
}

// process SAMTOOLS_INDEX {
//     //
// }

// process CALL_VARIANTS {
//     //
// }

// process FILTER_VARIANTS {
//     //
// }

// process INDEX_VARIANTS {
//     //
// }

// process GET_CONSENSUS {
//     //
// }

// process EXTRACT_CONSENSUS {
//     //
// }
