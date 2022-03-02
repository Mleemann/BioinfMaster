/*
* sraX module
*/

params.CONTAINER = "quay.io/biocontainers/srax:1.5--pl5262ha8f3691_1"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/srax:1.5--pl5262ha8f3691_1"
params.OUTPUT = "amrfinder_output"

process sraX_basic {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}")
    publishDir("resistance/sraX_basic", mode: 'copy')
    tag { fasta }
    //container params.CONTAINER
    conda '/scicore/home/egliadr/leeman0000/anaconda3/envs/sraX'

    input:
    tuple val (sample_id), path (fasta)

    output:
    path ("${sample_id}_sraX_basic.tsv"), emit: resistance
    path ("${sample_id}_sraX_basic/*"), emit: sraX_results

    script:
    """
    source /scicore/home/egliadr/leeman0000/.bashrc
    conda activate sraX
    mkdir input_srax
    cp ${fasta} input_srax
    sraX -i input_srax -o ${sample_id}_sraX_basic
    mv ${sample_id}_sraX_basic/Results/Summary_files/sraX_detected_ARGs.tsv ${sample_id}_sraX_basic.tsv
    """
}


process sraX_ext {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}")
    publishDir("resistance/sraX_ext", mode: 'copy')
    tag { fasta }
    //container params.CONTAINER
    conda '/scicore/home/egliadr/leeman0000/anaconda3/envs/sraX'

    input:
    tuple val (sample_id), path (fasta)

    output:
    path ("${sample_id}_sraX_ext.tsv"), emit: resistance
    path ("${sample_id}_sraX_ext/*"), emit: sraX_results

    script:
    """
    source /scicore/home/egliadr/leeman0000/.bashrc
    conda activate sraX
    mkdir input_srax
    cp ${fasta} input_srax
    sraX -db ext -i input_srax -o ${sample_id}_sraX_ext
    mv ${sample_id}_sraX_ext/Results/Summary_files/sraX_detected_ARGs.tsv ${sample_id}_sraX_ext.tsv
    """
}
