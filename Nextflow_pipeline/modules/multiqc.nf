/*
*  multiqc module
*/

params.CONTAINER = "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
params.OUTPUT = "multiqc_output"

process multiqc {
    // publishDir(params.OUTPUT, mode: 'copy')
    // publishDir("results/preprocessing/multiqc", mode: 'copy')
    container params.CONTAINER

    input:
    path (inputfiles)

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """
}
