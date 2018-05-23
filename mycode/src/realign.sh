#!/usr/bin/env bash

# Input files
OUT="EB_realign/$1"
tBAM="miniBAM/$2"  # 'miniBAM/9ee76dddd4e8194f4ed004397c494766.gencode.v19.U1U11.flank1k.nochr.bed.bam'
nBAM="miniBAM/$3"  # 'miniBAM/e6cdf28912d360a2b947f5911c79e214.gencode.v19.U1U11.flank1k.nochr.bed.bam'
BED="$HOME/mutational_analysis/files/gencode.v19.U1U11.nochr.bed"
# Software and scripts
MERGE_TAB="$HOME/mutational_analysis/scripts/merge_tab.R"  # merge tabs from bcftools query
BW2="$HOME/software/bowtie2-2.3.4.1-linux-x86_64/bowtie2"  # bowtie2
BAM2FQ="$HOME/mutational_analysis/scripts/bam2fq_for_realign.py"  # script for converting bam to fastq
BAMSTAT="$HOME/mutational_analysis/scripts/realign_bam_stat.py"

echo Run $1
#################
### PART 2
### Realign reads with MQ<30
#################
# Make FASTQ for realignments
python $BAM2FQ $tBAM $BED ${OUT}_tumour
python $BAM2FQ $nBAM $BED ${OUT}_normal
# Run Bowtie2
$BW2 --no-mixed --no-discordant -k 100 --very-sensitive -p 4 -x /db1/reference/genome -1 ${OUT}_tumour_1.fq -2 ${OUT}_tumour_2.fq -U ${OUT}_tumour_unpair.fq -S ${OUT}_tumour.sam
$BW2 --no-mixed --no-discordant -k 100 --very-sensitive -p 4 -x /db1/reference/genome -1 ${OUT}_normal_1.fq -2 ${OUT}_normal_2.fq -U ${OUT}_normal_unpair.fq -S ${OUT}_normal.sam
# sort and index
samtools sort ${OUT}_tumour.sam > ${OUT}_tumour.bam
samtools sort ${OUT}_normal.sam > ${OUT}_normal.bam
samtools index ${OUT}_tumour.bam
samtools index ${OUT}_normal.bam
# Get BAM STAT
python $BAMSTAT ${OUT}_tumour.bam ${OUT}_tumour.stat.csv
python $BAMSTAT ${OUT}_normal.bam ${OUT}_normal.stat.csv
