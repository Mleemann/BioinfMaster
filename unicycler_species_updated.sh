#!/bin/bash


#SBATCH --qos=1day
#SBATCH --mem=50G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=run3
#SBATCH --array=1-10


array=(mock 400119-16 600284-16 403902-15 500144-16 FR115 FR217 run210-49775 run306-49775 run369-49775 run410-49775)


source /scicore/home/egliadr/leeman0000/.bashrc
conda activate /scicore/home/egliadr/leeman0000/anaconda3/envs/Unicycler_updated

export _JAVA_OPTIONS="-Xmx10g"

export sample_id=${array["$SLURM_ARRAY_TASK_ID"]}

raw_R1="$PWD"/../reads/*/"$sample_id"_R1.fastq.gz
raw_R2="$PWD"/../reads/*/"$sample_id"_R2.fastq.gz

IFS='/' read -ra Path_reads <<< $(echo $raw_R1)
let "index = ${#Path_reads[@]} -2"
species=${Path_reads[$index]}

mkdir -p results/"$sample_id"/0_trimming

trimmomatic PE -threads "$SLURM_CPUS_PER_TASK" -phred33 $(echo "$raw_R1") $(echo "$raw_R2") results/"$sample_id"/0_trimming/r1.fastq.gz results/"$sample_id"/0_trimming/r1.not-paired.fastq.gz results/"$sample_id"/0_trimming/r2.fastq.gz results/"$sample_id"/0_trimming/r2.not-paired.fastq.gz ILLUMINACLIP:/scicore/home/egliadr/GROUP/Software/scripts/assembly_pipeline/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:12 MINLEN:100 2> results/"$sample_id"/0_trimming/quality_read_trimm_info

R1=results/"$sample_id"/0_trimming/r1.fastq.gz
R2=results/"$sample_id"/0_trimming/r2.fastq.gz

mkdir -p results/"$sample_id"/1_unicycler

unicycler -t "$SLURM_CPUS_PER_TASK" -1 "$R1" -2 "$R2" -o results/"$sample_id"/1_unicycler # optional long reads #  -l ../reads/"$sample_id".fastq.gz


#last pilon with  --mindepth 5 --minmq 10 option
mkdir -p results/"$sample_id"/1_unicycler/last_pilon
bwa index results/"$sample_id"/1_unicycler/assembly.fasta
bwa mem -t "$SLURM_CPUS_PER_TASK" results/"$sample_id"/1_unicycler/assembly.fasta "$R1" "$R2" > results/"$sample_id"/1_unicycler/last_pilon/alingnment.sam
samtools sort -@ "$SLURM_CPUS_PER_TASK" -T results/"$sample_id"/1_unicycler/last_pilon/temp_sort -o results/"$sample_id"/1_unicycler/last_pilon/alingnment.bam results/"$sample_id"/1_unicycler/last_pilon/alingnment.sam
samtools rmdup results/"$sample_id"/1_unicycler/last_pilon/alingnment.bam results/"$sample_id"/1_unicycler/last_pilon/alingnment.removed_duplicates.bam
samtools index results/"$sample_id"/1_unicycler/last_pilon/alingnment.removed_duplicates.bam
pilon --mindepth 5 --minmq 10 --threads "$SLURM_CPUS_PER_TASK" --genome results/"$sample_id"/1_unicycler/assembly.fasta --frags results/"$sample_id"/1_unicycler/last_pilon/alingnment.removed_duplicates.bam --changes --variant --outdir results/"$sample_id"/1_unicycler/last_pilon/pilon --output "$sample_id"


prokka --centre USB --compliant --addgenes --mincontiglen 200 --genus Genus --species species --prefix "$sample_id" --rfam --locustag "$sample_id" --strain "$sample_id" --outdir results/"$sample_id"/2_annotation --cpus "$SLURM_CPUS_PER_TASK" results/"$sample_id"/1_unicycler/last_pilon/pilon/"$sample_id".fasta

conda activate

conda activate /scicore/home/egliadr/leeman0000/anaconda3/envs/assembly_pipeline_updated
export PATH="/scicore/home/egliadr/GROUP/Software/scripts/raw_data_quality:$PATH"

mkdir -p results/"$sample_id"/3_quality

#---assembly stats---

quast --min-contig 0 -o results/"$sample_id"/3_quality/quast/ results/"$sample_id"/2_annotation/"$sample_id".fna

#---rMLST---

mkdir -p results/"$sample_id"/3_quality/rMLST/blast

counter=0
for gene in /scicore/home/egliadr/GROUP/Software/scripts/raw_data_quality/rMLST/database_20180323/*.fas
do
let counter=counter+1
blastn -num_threads "$SLURM_CPUS_PER_TASK" -db "$gene" -query results/"$sample_id"/2_annotation/"$sample_id".fna -max_target_seqs 100 -max_hsps 1   -outfmt "6 qseqid sseqid stitle qlen slen length pident nident mismatch gaps evalue bitscore" > results/"$sample_id"/3_quality/rMLST/blast/rMLST_blast_"$counter".tab

done
/scicore/home/egliadr/leeman0000/tools/python_scripts/call_rMLST_updated_P3.py /scicore/home/egliadr/GROUP/Software/scripts/raw_data_quality/rMLST/database_20180323/bigsdb-rMLST.txt results/"$sample_id"/3_quality/rMLST/blast/rMLST_blast*tab > results/"$sample_id"/3_quality/rMLST/"$sample_id"_rMLST.tab

#---Metaphlan3---

mkdir -p results/"$sample_id"/3_quality/Metaphlan3
#bwa mem -t "$SLURM_CPUS_PER_TASK" /scicore/home/egliadr/GROUP/Software/miniconda2/miniconda2/envs/raw_data_quality/bin/db_v20/mpa_v20_m200.fasta "$R1" "$R2" > results/"$sample_id"/3_quality/Metaphlan3/alignment.sam
bowtie2 --sam-no-hd --sam-no-sq --no-unal --very-sensitive -S results/"$sample_id"/3_quality/Metaphlan3/alignment.sam -x /scicore/home/egliadr/leeman0000/dbs/metaphlan/db/mpa_v30_CHOCOPhlAn_201901 -1 "$R1"  -2 "$R2"
metaphlan results/"$sample_id"/3_quality/Metaphlan3/alignment.sam --input_type sam --nproc "$SLURM_CPUS_PER_TASK" --index mpa_v30_CHOCOPhlAn_201901 --bowtie2db /scicore/home/egliadr/leeman0000/dbs/metaphlan/db/ > results/"$sample_id"/3_quality/Metaphlan3/profiled_metagenome.txt 2> results/"$sample_id"/3_quality/Metaphlan3/error.txt

#---remapping---

R1=results/"$sample_id"/0_trimming/r1.fastq.gz
R2=results/"$sample_id"/0_trimming/r2.fastq.gz

mkdir -p results/"$sample_id"/3_quality/remapping/mapping/
/scicore/home/egliadr/leeman0000/tools/python_scripts/make_one_contig_updated_P3.py results/"$sample_id"/2_annotation/"$sample_id".fna "$sample_id" > results/"$sample_id"/3_quality/remapping/"$sample_id"_concatenated_contigs.fna
bwa index results/"$sample_id"/3_quality/remapping/"$sample_id"_concatenated_contigs.fna results/"$sample_id"/3_quality/remapping/"$sample_id"_concatenated_contigs.fna
bwa mem -t "$SLURM_CPUS_PER_TASK" results/"$sample_id"/3_quality/remapping/"$sample_id"_concatenated_contigs.fna "$R1" "$R2" > results/"$sample_id"/3_quality/remapping/mapping/alingnment.sam
samtools sort -@ "$SLURM_CPUS_PER_TASK" -T results/"$sample_id"/3_quality/remapping/mapping/temp_sort -o results/"$sample_id"/3_quality/remapping/mapping/alingnment.bam results/"$sample_id"/3_quality/remapping/mapping/alingnment.sam
samtools rmdup results/"$sample_id"/3_quality/remapping/mapping/alingnment.bam results/"$sample_id"/3_quality/remapping/mapping/alingnment.removed_duplicates.bam
samtools index results/"$sample_id"/3_quality/remapping/mapping/alingnment.removed_duplicates.bam
/scicore/home/egliadr/leeman0000/tools/python_scripts/parse_sam_for_insertsize_updated_P3.py results/"$sample_id"/3_quality/remapping/mapping/alingnment.sam > results/"$sample_id"/3_quality/remapping/mapping/alingnment.removed_duplicates.bam.insertions.tab
/scicore/home/egliadr/GROUP/Software/scripts/raw_data_quality/insert_distribution.r results/"$sample_id"/3_quality/remapping/mapping/alingnment.removed_duplicates.bam.insertions.tab
mkdir -p results/"$sample_id"/3_quality/remapping/pilon
pilon --threads "$SLURM_CPUS_PER_TASK" --genome results/"$sample_id"/3_quality/remapping/"$sample_id"_concatenated_contigs.fna --frags results/"$sample_id"/3_quality/remapping/mapping/alingnment.removed_duplicates.bam --changes  --vcf --fix snps --outdir results/"$sample_id"/3_quality/remapping/pilon --output "$sample_id"
/scicore/home/egliadr/leeman0000/tools/python_scripts/make_coverage_pilon_corrected_updated_P3.py results/"$sample_id"/3_quality/remapping/pilon/"$sample_id".vcf > results/"$sample_id"/3_quality/remapping/pilon/coverage.tab
/scicore/home/egliadr/GROUP/Software/scripts/raw_data_quality/plot_coverage.r results/"$sample_id"/3_quality/remapping/pilon/coverage.tab

mv results/"$sample_id"/3_quality/remapping/mapping/alingnment.removed_duplicates.bam.insertions.tab results/"$sample_id"/3_quality/remapping/insertions.tab
mv results/"$sample_id"/3_quality/remapping/mapping/alingnment.removed_duplicates.bam.insertions.tab.png results/"$sample_id"/3_quality/remapping/insert_size.png
mv results/"$sample_id"/3_quality/remapping/pilon/coverage.tab results/"$sample_id"/3_quality/remapping/coverage.tab
mv results/"$sample_id"/3_quality/remapping/pilon/coverage.tab.png results/"$sample_id"/3_quality/remapping/coverage.png

#---find 16S species---
mkdir -p results/"$sample_id"/3_quality/16S
blastn -db /scicore/home/egliadr/GROUP/Software/scripts/assembly_pipeline/16S/db/16SMicrobial -num_threads "$SLURM_CPUS_PER_TASK" -max_target_seqs 1 -max_hsps 1 -query results/"$sample_id"/3_quality/remapping/"$sample_id"_concatenated_contigs.fna -out results/"$sample_id"/3_quality/16S/"$sample_id"_16S_blast.tab -outfmt "6 qseqid sseqid stitle qlen slen length pident nident mismatch gaps evalue bitscore"
echo -e NA'\t'NA'\t'NA'\t'NA'\t'NA'\t'NA'\t'NA >>  results/"$sample_id"/3_quality/16S/"$sample_id"_16S_blast.tab


#---search for virluence---
mkdir results/"$sample_id"/4_virulence/
#rgi -n "$SLURM_CPUS_PER_TASK" -i results/"$sample_id"/2_annotation/"$sample_id".faa -t protein -o results/"$sample_id"/4_virulence/"$sample_id"_resistence.tab -a BLAST
abricate --quiet --nopath --db ncbi  results/"$sample_id"/2_annotation/"$sample_id".fna > results/"$sample_id"/4_virulence/"$sample_id"_resistance.tab
abricate --quiet --nopath --db vfdb  results/"$sample_id"/2_annotation/"$sample_id".fna > results/"$sample_id"/4_virulence/"$sample_id"_virulence.tab

#---clean---

rm results/"$sample_id"/0_trimming/*.fastq.gz
rm results/"$sample_id"/3_quality/Metaphlan3/alignment.sam
rm results/"$sample_id"/3_quality/remapping/pilon/*
rm results/"$sample_id"/3_quality/remapping/mapping/*
rm results/"$sample_id"/3_quality/rMLST/blast/*
rm results/"$sample_id"/1_unicycler/last_pilon/*.bam
rm results/"$sample_id"/1_unicycler/last_pilon/*.sam

#---summary---

mkdir results/"$sample_id"/3_quality/summery
echo "$species" > results/"$sample_id"/3_quality/summery/9.txt
awk '/Input Read Pairs/{print $8}' results/"$sample_id"/0_trimming/quality_read_trimm_info | awk -F '%'  '{print $1}' | awk -F '('  '{print $2}' > results/"$sample_id"/3_quality/summery/1.txt
awk '/read_depth/{print $3}'  results/"$sample_id"/3_quality/remapping/coverage.tab | sort -n  | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }' > results/"$sample_id"/3_quality/summery/2.txt
grep -c alternative_base results/"$sample_id"/3_quality/remapping/coverage.tab > results/"$sample_id"/3_quality/summery/3.txt
grep -v Insert_size results/"$sample_id"/3_quality/remapping/insertions.tab | sort -n  | awk ' { a[i++]=$1; } END { print a[int(i/2)]; }' > results/"$sample_id"/3_quality/summery/4.txt
tail -n 1 results/"$sample_id"/3_quality/quast/transposed_report.tsv | awk '{print $14 "\t" $16  "\t" $18  "\t" $17 }' > results/"$sample_id"/3_quality/summery/5.txt
head -n 1 results/"$sample_id"/3_quality/16S/"$sample_id"_16S_blast.tab | awk -F "\t" '{print $3 "\t" $6  "\t" $7}'   > results/"$sample_id"/3_quality/summery/6.txt
grep "s__" results/"$sample_id"/3_quality/Metaphlan3/profiled_metagenome.txt  | grep -v "t__" | awk '{split($0,a,"|"); print a[7],"\t",$3}'| awk -F __ '{print $2}' | cut -f1,3 | head -n 1 > results/"$sample_id"/3_quality/summery/7.txt
cat results/"$sample_id"/3_quality/rMLST/"$sample_id"_rMLST.tab > results/"$sample_id"/3_quality/summery/8.txt

paste <(echo "$sample_id") results/"$sample_id"/3_quality/summery/*.txt > results/"$sample_id"/3_quality/summery/"$sample_id".tab

echo -e "Sample\tRead_quality\tRead_depth\tAlternative_bases\tInsert_size\tContig_count\tTotal_length\tN50\tGC_percent\t16S_species\tAlignment_length\tAlignment_idendity\tMetaPhlAn3_species\tMetaPhlAn3_purity\trMLST_best_species\trMLST_best_rST\tAlleles_missing\trMLST_2nd_best_rST\tAlleles_missing\trMLST_2nd_best_rST\tAlleles_missing\tinitial_species" > quality.tab
cat results/*/3_quality/summery/*.tab >> quality.tab
let sample_count=$(grep -c "" quality.tab)-1
echo "analysed_samples: $sample_count" >> quality.tab

sed 's/,/;/' quality.tab | sed 's/\t/,/g' > quality.csv



mkdir -p transfer_result/genomes/"$species"
mkdir -p transfer_result/genomes_one_contig/"$species"
mkdir -p transfer_result/reads/"$species"
mv quality.* transfer_result/

ln -s "$PWD"/results/"$sample_id"/2_annotation/"$sample_id".fna transfer_result/genomes/"$species"/
ln -s "$PWD"/results/"$sample_id"/3_quality/remapping/"$sample_id"_concatenated_contigs.fna transfer_result/genomes_one_contig/"$species"/
ln -s $(echo "$raw_R1") transfer_result/reads/"$species"/
ln -s $(echo "$raw_R2") transfer_result/reads/"$species"/
