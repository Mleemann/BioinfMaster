/*
* RGI module
*/

params.CONTAINER = "quay.io/biocontainers/rgi:5.2.0--pyhdfd78af_0"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/rgi:5.2.0--pyhdfd78af_0"
params.OUTPUT = "rgi_output"

process rgi {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}")
    publishDir("resistance/rgi", mode: 'copy')
    tag { fasta }
    //container params.CONTAINER
    conda '/scicore/home/egliadr/leeman0000/anaconda3/envs/rgi520'

    input:
    tuple val (sample_id), path (fasta)
    path (db)

    output:
    path ("${sample_id}_rgi_resistance.tab.txt"), emit: full_table
    path ("${sample_id}_rgi.tab"), emit: resistance

    script:
    """
    source /scicore/home/egliadr/leeman0000/.bashrc
    conda activate rgi520
    rgi main -i ${fasta} -t contig -o ${sample_id}_rgi_resistance.tab -a BLAST --clean --local
    cat ${sample_id}_rgi_resistance.tab.txt | cut -f2,9,10,11,12,14,15,16 > ${sample_id}_rgi.tab
    """
}
