#!/bin/bash


#SBATCH --qos=30min
#SBATCH --mem=5G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=get_AMR_data
#SBATCH --time=00:10:00

species="KP"

while read -r file; do grep -v "#" */resistance/abricate/"$file"_abricate_ncbi.tab >> "$species"_abricate_ncbi.tab; done < samples_"$species".txt
while read -r file; do grep -v "#" */resistance/abricate/"$file"_abricate_card.tab >> "$species"_abricate_card.tab; done < samples_"$species".txt
while read -r file; do grep -v "#" */resistance/abricate/"$file"_abricate_megares.tab >> "$species"_abricate_megares.tab; done < samples_"$species".txt
while read -r file; do grep -v "#" */resistance/abricate/"$file"_abricate_resfinder.tab >> "$species"_abricate_resfinder.tab; done < samples_"$species".txt
while read -r file; do grep -v "#" */resistance/abricate/"$file"_abricate_argannot.tab >> "$species"_abricate_argannot.tab; done < samples_"$species".txt


while read -r file; do awk FNR-1 */resistance/amrfinder_nuc/"$file"_amrfinder_nuc.tsv >> "$species"_amrfinder_nuc.tab; done < samples_"$species".txt
while read -r file; do awk FNR-1 */resistance/amrfinder_prot/"$file"_amrfinder_prot.tsv >> "$species"_amrfinder_prot.tab; done < samples_"$species".txt


while read -r file; do awk FNR-1 */resistance/deeparg_LS/"$file"_deeparg_LS.tsv.mapping*.ARG >> "$species"_deeparg_LS.tab; done < samples_"$species".txt
grep -v .potential. "$species"_deeparg_LS.tab > "$species"_deeparg_LS_noPot.tab

# use grep to include from which file the entries are
while read -r file; do grep -v "#" */resistance/deeparg_SR/"$file"_deeparg_SR.tsv.clean.deeparg.mapping*.ARG >> "$species"_deeparg_SR.tab; done < samples_"$species".txt
grep -v .potential. "$species"_deeparg_SR.tab > "$species"_deeparg_SR_noPot.tab

while read -r file; do awk FNR-1 */resistance/resfinder_assembly/"$file"_resfinder_assembly/ResFinder_results_tab.txt >> "$species"_resfinder_assembly.tab; done < samples_"$species".txt

# to add sample_id to each entry
while read -r file; do sed -i "s/$/\t$file/" */resistance/resfinder_reads/"$file"_resfinder_reads/ResFinder_results_tab.txt; done < samples_"$species".txt
while read -r file; do awk FNR-1 */resistance/resfinder_reads/"$file"_resfinder_reads/ResFinder_results_tab.txt >> "$species"_resfinder_reads.tab; done < samples_"$species".txt

while read -r file; do cat */resistance/rgi/"$file"_rgi_resistance.tab.txt | cut -f2,9,10,11,12,14,15,16,17 >> "$species"_rgi.tab; done < samples_"$species".txt

while read -r file; do awk FNR-1 */resistance/sraX_basic/"$file"_sraX_basic/Results/Summary_files/sraX_gene_coordinates.tsv >> "$species"_srax_basic.tab; done < samples_"$species".txt
while read -r file; do awk FNR-1 */resistance/sraX_ext/"$file"_sraX_ext/Results/Summary_files/sraX_gene_coordinates.tsv >> "$species"_srax_ext.tab; done < samples_"$species".txt
