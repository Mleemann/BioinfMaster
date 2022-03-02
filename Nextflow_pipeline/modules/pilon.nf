/*
*  pilon module
*/

params.CONTAINER = "quay.io/biocontainers/pilon:1.24--hdfd78af_0"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/pilon:1.24--hdfd78af_0"
params.OUTPUT = ""

process pilon {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}/polished_assembly", mode: 'copy')
    tag { "${sample_id}" }
    container params.CONTAINER

    input:
    tuple val (sample_id), path(bam), path(assembly)

    output:
    tuple val (sample_id), path("${sample_id}.fasta"), emit: assembly
    tuple val (sample_id), path("${sample_id}.vcf"), emit: vcf

    script:
    """
    export _JAVA_OPTIONS="-Xmx10g"
    pilon --mindepth 5 --minmq 10 --threads ${task.cpus} \
    --genome ${assembly} --frags ${bam} --changes --variant \
    --outdir ./ --output ${sample_id}
    """
}


process pilon_remapping {
    // publishDir(params.OUTPUT, mode: 'copy')
    tag { "${sample_id}" }
    container params.CONTAINER

    input:
    tuple val (sample_id), path(bam), path(assembly)

    output:
    tuple val (sample_id), path("${sample_id}.fasta"), emit: assembly
    tuple val (sample_id), path("${sample_id}.vcf"), emit: vcf

    script:
    """
    export _JAVA_OPTIONS="-Xmx10g"
    pilon --threads ${task.cpus} --genome ${assembly} --frags ${bam} \
    --changes --vcf --fix snps --outdir ./ --output ${sample_id}
    """
}
