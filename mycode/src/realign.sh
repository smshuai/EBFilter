#!/usr/bin/env bash

# Input arguments
miniBAM=$1  # 'miniBAM/e6cdf28912d360a2b947f5911c79e214.gencode.v19.U1U11.flank1k.nochr.bed.bam'
donorID=$2  # donor/patient identifier
OUT=$3  # output prefix
# References
BW2IX="/db1/reference/genome"  # bowtie2 index
U1BED=""  # U1 bed with all variants and pseudogenes
U1CORE=""  # Core U1 genes (7 U1 + 1 U11) 
# Software and references
BW2="$HOME/software/bowtie2-2.3.4.1-linux-x86_64/bowtie2"  # bowtie2
src_dir=""  # Directory for scripts if not in PATH
BAM2FQ="$src_dir/bam2fq_for_realign.py"  # script for converting bam to fastq
POSTPROCESS="$src_dir/postprocess_realign.py"  # script for postprocessing
CALCDP="$src_dir/calc_weighted_depth.py"  # script used to calcuate weighted depth

#################
### PART 1
### Realign reads with MQ<30
#################
# Make FASTQ for realignments
python $BAM2FQ $miniBAM $U1BED $OUT
# Run Bowtie2
$BW2 --no-mixed --no-discordant -k 100 --very-sensitive -p 4 -x $BW2IX -1 ${OUT}_1.fq -2 ${OUT}_2.fq -U ${OUT}_unpair.fq -S ${OUT}_realign.sam
# sort and index
samtools sort ${OUT}_realign.sam > ${OUT}_realign.bam
samtools index ${OUT}_realign.bam
#################
### PART 2
### Post-process re-aligned SAM
#################
python $POSTPROCESS ${OUT}_realign.bam $miniBAM $OUT $U1BED
#################
### PART 3
### Make weighted depth table
#################
python $CALCDP ${OUT}_processed.bam ${OUT}_mstat.tsv ${OUT}_umap.txt $U1CORE $donorID ${OUT}_weighted_dp.tsv
## Clean up
rm ${OUT}_1.fq ${OUT}_2.fq ${OUT}_unpair.fq ${OUT}_realign.sam