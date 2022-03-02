/*
*  bwa index module
*/

params.CONTAINER = "quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/bwa:0.7.17--h5bf99c6_8"
params.OUTPUT = ""

process bwaIndex {
    tag { sample_id }
    container params.CONTAINER

    input:
    tuple val (sample_id), path(assembly)

    output:
    tuple val (sample_id), path("${assembly}.*"), emit: index

    script:
    """
    bwa index ${assembly}
    """
}
