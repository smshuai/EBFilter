#!/usr/bin/env bash

# Input files
miniBAM=$1  # 'miniBAM/e6cdf28912d360a2b947f5911c79e214.gencode.v19.U1U11.flank1k.nochr.bed.bam'
BED="$HOME/mutational_analysis/files/gencode.v19.U1U11.nochr.bed"
OUT=$2
# Software and scripts
BW2="$HOME/software/bowtie2-2.3.4.1-linux-x86_64/bowtie2"  # bowtie2
BAM2FQ="$HOME/mutational_analysis/scripts/bam2fq_for_realign.py"  # script for converting bam to fastq

echo Run $1
#################
### PART 2
### Realign reads with MQ<30
#################
# Make FASTQ for realignments
python $BAM2FQ $miniBAM $BED $OUT
# Run Bowtie2
$BW2 --no-mixed --no-discordant -k 100 --very-sensitive -p 4 -x /db1/reference/genome -1 ${OUT}_1.fq -2 ${OUT}_2.fq -U ${OUT}_unpair.fq -S ${OUT}_realign.sam
# sort and index
samtools sort ${OUT}_realign.sam > ${OUT}_realign.bam
samtools index ${OUT}_realign.bam