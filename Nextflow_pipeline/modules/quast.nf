/*
*  quast module
*/

params.CONTAINER = "quay.io/biocontainers/quast:5.0.2--py37pl5262h190e900_4"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl5262h190e900_4"
params.OUTPUT = "quast_output"

process quast {
    // publishDir(params.OUTPUT, mode: 'copy')
    tag { fasta }
    container params.CONTAINER

    input:
    tuple val (sample_id), path (fasta)

    output:
    tuple val (sample_id), path ("${sample_id}/transposed_report.tsv"), emit: tsv
    path ("${sample_id}")

    script:
    """
    quast --min-contig 0 -o ${sample_id} ${fasta}
    """
}
