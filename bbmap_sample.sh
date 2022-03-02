#!/bin/bash


#SBATCH --qos=30min
#SBATCH --mem=3G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=Downsampling
#SBATCH --time=00:30:00

module load BBMap/38.91-foss-2018b

for reads in 217740 206280 194820 183360
do
reformat.sh in1="../../reads/staaur/129011-20_R1.fastq.gz" in2="../../reads/staaur/129011-20_R2.fastq.gz" out1="staaur/129011-20_${reads}x_R1.fastq.gz" out2="staaur/129011-20_${reads}x_R2.fastq.gz" samplereadstarget=$reads sampleseed=30
done

for reads in 160681 149204 137727 126250 114772
do
reformat.sh in1="../../reads/staaur/351682-21_R1.fastq.gz" in2="../../reads/staaur/351682-21_R2.fastq.gz" out1="staaur/351682-21_${reads}x_R1.fastq.gz" out2="staaur/351682-21_${reads}x_R2.fastq.gz" samplereadstarget=$reads sampleseed=30
done

for reads in 224881 213045 201210 189374
do
reformat.sh in1="../../reads/staaur/360552-21_R1.fastq.gz" in2="../../reads/staaur/360552-21_R2.fastq.gz" out1="staaur/360552-21_${reads}x_R1.fastq.gz" out2="staaur/360552-21_${reads}x_R2.fastq.gz" samplereadstarget=$reads sampleseed=30
done

for reads in 157074 145854 134635 123415 112195
do
reformat.sh in1="../../reads/mrsa/402892-20_R1.fastq.gz" in2="../../reads/mrsa/402892-20_R2.fastq.gz" out1="mrsa/402892-20_${reads}x_R1.fastq.gz" out2="mrsa/402892-20_${reads}x_R2.fastq.gz" samplereadstarget=$reads sampleseed=30
done

for reads in 216395 205006 193617 182228
do
reformat.sh in1="../../reads/mrsa/802440-21_R1.fastq.gz" in2="../../reads/mrsa/802440-21_R2.fastq.gz" out1="mrsa/802440-21_${reads}x_R1.fastq.gz" out2="mrsa/802440-21_${reads}x_R2.fastq.gz" samplereadstarget=$reads sampleseed=30
done
