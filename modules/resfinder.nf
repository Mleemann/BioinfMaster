/*
* resfinder module
*/

params.CONTAINER = "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
params.OUTPUT = "resfinder_output"

process resfinder_fasta {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}")
    publishDir("resistance/resfinder_assembly", mode: 'copy')
    tag { sample_id }
    //container params.CONTAINER
    conda '/scicore/home/egliadr/leeman0000/anaconda3/envs/blast_python'

    input:
    tuple val (sample_id), path (fasta)
    path (db)
    val (species)

    output:
    path ("${sample_id}_resfinder_assembly/*"), emit: resistance

    script:
    """
    source /scicore/home/egliadr/leeman0000/.bashrc
    conda activate blast_python
    python3 /scicore/home/egliadr/leeman0000/tools/resfinder/run_resfinder.py \
    -o ${sample_id}_resfinder_assembly -s "${species}" \
    -l 0.6 -t 0.8 --acquired -db_res ${db} -ifa ${fasta}
    """
}

process resfinder_reads {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}")
    publishDir("resistance/resfinder_reads", mode: 'copy')
    tag { sample_id }
    //container params.CONTAINER
    conda '/scicore/home/egliadr/leeman0000/anaconda3/envs/blast_python'

    input:
    tuple val (sample_id), path (fastq_r1), path (fastq_r2)
    path (db)
    val (species)

    output:
    path ("${sample_id}_resfinder_reads/*"), emit: resistance

    script:
    """
    source /scicore/home/egliadr/leeman0000/.bashrc
    conda activate blast_python
    python3 /scicore/home/egliadr/leeman0000/tools/resfinder/run_resfinder.py \
    -o ${sample_id}_resfinder_reads -s "${species}" \
    -l 0.6 -t 0.8 --acquired -db_res ${db} -ifq ${fastq_r1} ${fastq_r2}
    """
}
