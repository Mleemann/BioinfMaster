/*
*  samtools module
*/

params.CONTAINER = "quay.io/biocontainers/samtools:1.13--h8c37831_0"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/samtools:1.13--h8c37831_0"
params.OUTPUT = ""

process samtools {
    tag { sample_id }
    container params.CONTAINER

    input:
    tuple val (sample_id), path(sam)

    output:
    tuple val (sample_id), path("*.removed_duplicates.bam"), emit: bam

    script:
    """
    samtools sort -@ ${task.cpus} -T ${sample_id} -o ${sample_id}_alingnment.bam ${sam}
    samtools rmdup ${sample_id}_alingnment.bam ${sample_id}_alingnment.removed_duplicates.bam
    samtools index ${sample_id}_alingnment.removed_duplicates.bam
    """
}
