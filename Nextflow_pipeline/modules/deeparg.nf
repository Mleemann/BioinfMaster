/*
* deeparg module
*/

params.CONTAINER = "docker pull gaarangoa/deeparg:latest"
params.OUTPUT = "deeparg_output"

process deeparg_LS {
    publishDir("results/${sample_id}/deeparg_LS")
    publishDir("resistance/deeparg_LS", mode: 'copy')
    tag { sample_id }
    //container params.CONTAINER
    conda '/scicore/home/egliadr/leeman0000/anaconda3/envs/deeparg_env'

    input:
    tuple val (sample_id), path (fasta)
    path (db)

    output:
    path ("${sample_id}_deeparg_LS.*.ARG"), emit: resistance

    script:
    """
    source /scicore/home/egliadr/leeman0000/.bashrc
    conda activate deeparg_env
    deeparg predict --model LS --type nucl --input ${fasta} -d ${db} \\
            --arg-alignment-identity 80 --out ${sample_id}_deeparg_LS.tsv
    """
}

process deeparg_SR {
    publishDir("results/${sample_id}/deeparg_SR")
    publishDir("resistance/deeparg_SR", mode: 'copy')
    tag { sample_id }
    //container params.CONTAINER
    conda '/scicore/home/egliadr/leeman0000/anaconda3/envs/deeparg_env'

    input:
    tuple val (sample_id), path (fastq_r1), path (fastq_r2)
    path (db)

    output:
    path ("${sample_id}_deeparg_SR.tsv*"), emit: resistance

    script:
    """
    source /scicore/home/egliadr/leeman0000/.bashrc
    conda activate deeparg_env
    deeparg short_reads_pipeline --forward_pe_file ${fastq_r1} --reverse_pe_file ${fastq_r2} \\
          --output_file ${sample_id}_deeparg_SR.tsv -d ${db} --bowtie_16s_identity 100
    """
}
