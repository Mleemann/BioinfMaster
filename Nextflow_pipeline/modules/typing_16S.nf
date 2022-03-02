/*
*  16S module
*/

params.CONTAINER = "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
//params.CONTAINER = "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
params.OUTPUT = "typing16s_output"

process typing_16S {
    // publishDir(params.OUTPUT, mode: 'copy')
    tag { sample_id }
    container params.CONTAINER

    input:
    //tuple val (sample_id), path (one_contig)
    tuple val (sample_id), path (one_contig)
    val (db)

    output:
    tuple val (sample_id), path ("${sample_id}_16S_blast.tab"), emit: blast_tab

    script:
    """
    #!/bin/bash
    DB=`find -L ${db} -name "*.naa" | sed 's/.naa//'`
    blastn -db \$DB  -num_threads ${task.cpus} -max_target_seqs 1 -max_hsps 1 \
           -query ${one_contig} -out ${sample_id}_16S_blast.tab \
           -outfmt "6 qseqid sseqid stitle qlen slen length pident nident mismatch gaps evalue bitscore"
    echo -e NA'\t'NA'\t'NA'\t'NA'\t'NA'\t'NA'\t'NA >> ${sample_id}_16S_blast.tab
    """
}
