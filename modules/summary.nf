/*
*  summary module
*/

params.OUTPUT = "summary"

process summary_sample {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("results/${sample_id}/summary", mode: 'copy')
    tag { sample_id }

    input:
    tuple val (sample_id), path (trim_log), path (coverage), path (insertions), path (quast_tsv), path (blast16S), path (metaphlan), path (rmlst)

    output:
    path ("*.txt"), emit: summary_files
    path ("${sample_id}.tab"), emit: sample_quality

    script:
    """
    #!/bin/bash
    awk '/Input Read Pairs/{print \$8}' ${trim_log} | awk -F '%'  '{print \$1}' | awk -F '('  '{print \$2}' > ${sample_id}_1.txt
    awk '/read_depth/{print \$3}' ${coverage} | sort -n  | awk ' { a[i++]=\$1; } END { print a[int(i/2)]; }' > ${sample_id}_2.txt
    grep -c alternative_base ${coverage} > ${sample_id}_3.txt
    grep -v Insert_size ${insertions} | sort -n  | awk ' { a[i++]=\$1; } END { print a[int(i/2)]; }' > ${sample_id}_4.txt
    tail -n 1 ${quast_tsv} | awk '{print \$14 "\\t" \$16  "\\t" \$18  "\\t" \$17 }' > ${sample_id}_5.txt
    head -n 1 ${blast16S} | awk -F "\\t" '{print \$3 "\\t" \$6  "\\t" \$7}'   > ${sample_id}_6.txt
    grep "s__" ${metaphlan}  | grep -v "t__" | awk '{split(\$0,a,"|"); print a[7],"\\t",\$3}'| awk -F __ '{print \$2}' | cut -f1,3 | head -n 1 > ${sample_id}_7.txt
    cat ${rmlst} > ${sample_id}_8.txt
    cat <(echo ${sample_id}) ${sample_id}_?.txt | tr "\n" "\t" > ${sample_id}.tab
    """
}


process merge_summaries {
    // publishDir(params.OUTPUT, mode: 'copy')
    publishDir("summary", mode: 'copy')

    input:
    path (sample_quality)

    output:
    path ("quality.*"), emit: quality

    script:
    """
    #!/bin/bash
    echo -e "Sample\\tRead_quality\\tRead_depth\\tAlternative_bases\\tInsert_size\\tContig_count\\tTotal_length\\tN50\\tGC_percent\t16S_species\tAlignment_length\\tAlignment_idendity\\tMetaPhlAn3_species\\tMetaPhlAn3_purity\\trMLST_best_species\\trMLST_best_rST\\tAlleles_missing\\trMLST_2nd_best_rST\\tAlleles_missing\\trMLST_2nd_best_rST\\tAlleles_missing\\tinitial_species" > quality.tab
    for sample in ${sample_quality}; do cat \$sample >> quality.tab; printf "\n" >> quality.tab; done
    let sample_count=\$(grep -c "" quality.tab)-1
    echo "analysed_samples: \$sample_count" >> quality.tab
    sed 's/,/;/' quality.tab | sed 's/\t/,/g' > quality.csv
    """
}
