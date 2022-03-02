/*
*  abricate module
*/

params.CONTAINER = "quay.io/biocontainers/abricate:1.0.1--ha8f3691_1"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/abricate:1.0.1--ha8f3691_1"
params.OUTPUT = "abricate_output"

process abricate {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}")
    publishDir("resistance/abricate", mode: 'copy')
    tag { fasta }
    container params.CONTAINER

    input:
    tuple val (sample_id), path (fasta)
    val (db)

    output:
    path ("${sample_id}_abricate_*.tab"), emit: resistance

    script:
    """
    #!/bin/bash
    echo ${db}
    for database in ${db}; do abricate --quiet --nopath --db "\$database" ${fasta} > ${sample_id}_abricate_"\$database".tab; done
    """
}
