/*
*  bwa-mem module
*/

params.CONTAINER = "quay.io/biocontainers/bwa:0.7.17--h5bf99c6_8"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/bwa:0.7.17--h5bf99c6_8"
params.OUTPUT = ""

process bwaAlign{
    tag { sample_id }
    container params.CONTAINER

    input:
    tuple val (sample_id), path (fastq_r1), path (fastq_r2), path (index)

    output:
    tuple val (sample_id), path ("${sample_id}_alignment.sam"), emit: sam

    script:
    def indexname = index[0].baseName

    """
    bwa mem -t ${task.cpus} ${indexname} ${fastq_r1} ${fastq_r2} > ${sample_id}_alignment.sam
    """
}
