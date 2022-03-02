/*
*  unicycler module
*/

params.CONTAINER = "quay.io/biocontainers/unicycler:0.4.8--py37h13b99d1_3"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/unicycler:0.4.8--py37h13b99d1_3"
params.OUTPUT = "unicycler_output"

process unicycler {
    // publishDir(params.OUTPUT, mode: 'copy')
    // publishDir("results/${sample_id}/1_unicycler", mode: 'copy')
    tag { sample_id }
    container params.CONTAINER

    input:
    tuple val (sample_id), path (fastq_r1), path (fastq_r2)

    output:
    tuple val (sample_id), path ("${sample_id}_assembly.fasta"), emit: assembly
    tuple val (sample_id), path ("${sample_id}_unicyler.log"), emit: log

    script:
    """
    unicycler -t ${task.cpus} \
    -1 ${fastq_r1} -2 ${fastq_r2} \
    -o ./
    mv assembly.fasta ${sample_id}_assembly.fasta
    mv unicycler.log ${sample_id}_unicyler.log
    """
}
