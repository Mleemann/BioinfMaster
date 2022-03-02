/*
*  python functions
*/

params.OUTPUT = ""

process make_one_contig {
    // publishDir(params.OUTPUT, mode: 'copy')
    // publishDir("results/${sample_id}/3_quality/remapping", mode: 'copy')
    tag { fasta }

    input:
    tuple val (sample_id), path (fasta)

    output:
    tuple val (sample_id), path ("${sample_id}_concatenated_contigs.fna"), emit: one_contig

    script:
    """
    $baseDir/python_scripts/make_one_contig_updated_P3.py ${fasta} ${sample_id} \
    > ${sample_id}_concatenated_contigs.fna
    """
}


process parse_sam_for_insertsize {
    // publishDir(params.OUTPUT, mode: 'copy')
    // publishDir("results/${sample_id}/3_quality/remapping", mode: 'copy')
    tag { sam }

    input:
    tuple val (sample_id), path (sam)

    output:
    tuple val (sample_id), path ("${sample_id}.insertions.tab"), emit: insertions_tab

    script:
    """
    $baseDir/python_scripts/parse_sam_for_insertsize_updated_P3.py ${sam} \
    > ${sample_id}.insertions.tab
    """
}


process coverage_pilon_corrected {
    // publishDir(params.OUTPUT, mode: 'copy')
    // publishDir("results/${sample_id}/3_quality/remapping", mode: 'copy')
    tag { vcf }

    input:
    tuple val (sample_id), path (vcf)

    output:
    tuple val (sample_id), path ("${sample_id}_coverage.tab"), emit: coverage_tab

    script:
    """
    $baseDir/python_scripts/make_coverage_pilon_corrected_updated_P3.py ${vcf} \
    > ${sample_id}_coverage.tab
    """
}
