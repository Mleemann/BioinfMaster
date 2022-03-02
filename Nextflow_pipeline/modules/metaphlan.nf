/*
*  metaphlan module
*/

params.CONTAINER = "quay.io/biocontainers/metaphlan:3.0.13--pyhb7b1952_0"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/metaphlan:3.0.13--pyhb7b1952_0"
params.OUTPUT = "metaphlan_output"


process metaphlan3 {
    // publishDir(params.OUTPUT, mode: 'copy')
    tag { sample_id }
    container params.CONTAINER

    input:
    tuple val (sample_id), path (fastq_r1), path (fastq_r2)

    output:
    tuple val (sample_id), path ("${sample_id}_profiled_metagenome.txt"), emit: profile
    tuple val (sample_id), path ("${sample_id}.error.txt"), emit: error

    script:
    """
    bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive  \
            -S ${sample_id}_alignment.sam \
            -x /scicore/home/egliadr/leeman0000/dbs/metaphlan/db/mpa_v30_CHOCOPhlAn_201901 \
            -1 ${fastq_r1} -2 ${fastq_r2}
    metaphlan ${sample_id}_alignment.sam --input_type sam --nproc ${task.cpus} \
              --index mpa_v30_CHOCOPhlAn_201901 \
              --bowtie2db /scicore/home/egliadr/leeman0000/dbs/metaphlan/db/ \
              > ${sample_id}_profiled_metagenome.txt \
              2> ${sample_id}.error.txt
    """
}
