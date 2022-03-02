/*
*  prokka module
*/

 // "bioconda::prokka=1.14.6"
params.CONTAINER = "quay.io/biocontainers/prokka:1.14.6--pl5262hdfd78af_1"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl5262hdfd78af_1"
params.OUTPUT = "prokka_output"

process prokka {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}/annotation", mode: 'copy')
    tag { sample_id }
    container params.CONTAINER


    input:
    tuple val (sample_id), path (fasta)

    output:
    tuple val(sample_id), path ("${sample_id}/${sample_id}.fna"), emit: fna
    tuple val(sample_id), path ("${sample_id}/${sample_id}.faa"), emit: faa

    script:
    """
    prokka \\
      --centre USB --compliant --addgenes --mincontiglen 200 --genus Genus \\
      --species species --prefix ${sample_id} --rfam --locustag ${sample_id} \\
      --strain ${sample_id} --cpus ${task.cpus} ${fasta}
    """
}
