/*
*  trimmomatic module
*/

params.CONTAINER = "quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/trimmomatic:0.39--hdfd78af_2"
params.OUTPUT = "trimmomatic_output"


process trimmomaticPE {
    // publishDir(params.OUTPUT, mode: 'copy')
    // publishDir("results/${sample_id}/0_trimming", mode: 'copy')
    tag { sample_id }
    container params.CONTAINER

    input:
    tuple val (sample_id), path (fastq)
    val (illuminaclip)

    output:
    tuple val (sample_id), path ("${sample_id}_r1.fastq.gz"), path ("${sample_id}_r2.fastq.gz"), emit: trimmed_reads
    tuple val (sample_id), path ("${sample_id}.quality_read_trimm_info"), emit: trim_log

    script:
    """
    trimmomatic PE \
    -threads ${task.cpus} -phred33 \
    ${fastq} \
    ${sample_id}_r1.fastq.gz \
    ${sample_id}_r1.not-paired.fastq.gz \
    ${sample_id}_r2.fastq.gz \
    ${sample_id}_r2.not-paired.fastq.gz \
    ILLUMINACLIP:${illuminaclip}:2:30:10 SLIDINGWINDOW:4:12 MINLEN:100 \
    2> ${sample_id}.quality_read_trimm_info
    """
}
