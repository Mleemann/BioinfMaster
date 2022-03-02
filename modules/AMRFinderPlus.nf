/*
* AMRFinderPlus module
*/

params.CONTAINER = "quay.io/biocontainers/ncbi-amrfinderplus:3.10.18--h17dc2d4_0"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.10.18--h17dc2d4_0"
params.OUTPUT = "amrfinder_output"

process amrfinder_nuc {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}")
    publishDir("resistance/amrfinder_nuc", mode: 'copy')
    tag { fasta }
    //container params.CONTAINER
    conda '/scicore/home/egliadr/leeman0000/anaconda3/envs/amrfinder'

    input:
    tuple val (sample_id), path (fasta)

    output:
    path ("${sample_id}_amrfinder_nuc.tsv"), emit: resistance

    script:
    """
    source /scicore/home/egliadr/leeman0000/.bashrc
    conda activate amrfinder
    amrfinder -n ${fasta} -o ${sample_id}_amrfinder_nuc.tsv
    """
}


process amrfinder_prot {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}")
    publishDir("resistance/amrfinder_prot", mode: 'copy')
    tag { fasta }
    //container params.CONTAINER
    conda '/scicore/home/egliadr/leeman0000/anaconda3/envs/amrfinder'

    input:
    tuple val (sample_id), path (faa)

    output:
    path ("${sample_id}_amrfinder_prot.tsv"), emit: resistance

    script:
    """
    source /scicore/home/egliadr/leeman0000/.bashrc
    conda activate amrfinder
    amrfinder -p ${faa} -o ${sample_id}_amrfinder_prot.tsv
    """
}
