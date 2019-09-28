#!/bin/bash

# exit when any command fails
set -e

# script for aligning the fastq files to the oligos
# needs the paths below to minimap2, samtools and the directory containing the oligos and fastq files
# we used samtools 1.9 (https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2)
# and minimap2 version 2.13 (https://github.com/lh3/minimap2/releases/tag/v2.13)
MINIMAP2='/raid/shubham/minimap2/minimap2'
SAMTOOLS='/raid/shubham/samtools-1.9/samtools'
DATA_PATH='../../LDPC_DNA_storage_data_genapsys/'

## create .bed files from .fa files (will be used later for separating reads from different experiments)
#for i in {1..3}; do
#    # first create fai file
#    $SAMTOOLS faidx $DATA_PATH/oligo_files/oligos_$i.fa
#    awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' $DATA_PATH/oligo_files/oligos_$i.fa.fai > $DATA_PATH/oligo_files/oligos_$i.bed
#done
#
## align fastq
#mkdir -p $DATA_PATH/aligned
#for i in {1..3}; do
#    $MINIMAP2 -ax sr $DATA_PATH/oligo_files/oligos_$i.fa $DATA_PATH/fastq/exp_$i.fastq > $DATA_PATH/aligned/exp_$i.sam
#done
#
## now separate aligned sam file 
#for i in {1..3}; do
#    $SAMTOOLS view -h -L $DATA_PATH/oligo_files/oligos_$i.bed $DATA_PATH/aligned/exp_$i.sam > $DATA_PATH/aligned/exp_aligned_$i.sam
#done
#
## generate FASTQ and read files for each experiment (to be used for decoding)
#for i in {1..3}; do
#    $SAMTOOLS fastq -0 $DATA_PATH/aligned/exp_aligned_$i.fastq --reference $DATA_PATH/oligo_files/oligos_$i.fa $DATA_PATH/aligned/exp_aligned_$i.sam
#    sed -n '2~4p' $DATA_PATH/aligned/exp_aligned_$i.fastq > $DATA_PATH/aligned/exp_aligned_$i.reads
#done
#
## convert to BAM and sort (for collecting stats)
#for i in {1..3}; do
#    $SAMTOOLS view -b -o $DATA_PATH/aligned/exp_aligned_$i.bam $DATA_PATH/aligned/exp_aligned_$i.sam
#    $SAMTOOLS sort -O BAM -o $DATA_PATH/aligned/exp_aligned_sorted_$i.bam --reference $DATA_PATH/oligo_files/oligos_$i.fa $DATA_PATH/aligned/exp_aligned_$i.bam
#done

# run samtools stats for computing error rates 
mkdir -p $DATA_PATH/stats
for i in {2..3}; do
    $SAMTOOLS stats -r $DATA_PATH/oligo_files/oligos_$i.fa $DATA_PATH/aligned/exp_aligned_sorted_$i.bam > $DATA_PATH/stats/exp_aligned_sorted_$i.stats
    grep ^MPC $DATA_PATH/stats/exp_aligned_sorted_$i.stats | cut -f 2- > $DATA_PATH/stats/exp_aligned_sorted_$i.substitution_stats
    grep ^IC $DATA_PATH/stats/exp_aligned_sorted_$i.stats | cut -f 2,4 > $DATA_PATH/stats/exp_aligned_sorted_$i.insertion_stats 
    grep ^IC $DATA_PATH/stats/exp_aligned_sorted_$i.stats | cut -f 2,6 > $DATA_PATH/stats/exp_aligned_sorted_$i.deletion_stats 
    num_mapped=$(grep "reads mapped:" $DATA_PATH/stats/exp_aligned_sorted_$i.stats | cut -f 3)
    python3 ../util/compile_plot_stats.py $DATA_PATH/stats/exp_aligned_sorted_$i $num_mapped
done

# compute coverage stats (subsample to 5x and compute the coefficient of variation)
subsampling_coverage=5
for i in {1..3}; do
    # first extract the ref name from sam
    awk '$1 !~ "^@"' < $DATA_PATH/aligned/exp_aligned_$i.sam | cut -f 3 > $DATA_PATH/aligned/exp_aligned_$i.ref
    num_oligos=$(wc -l $DATA_PATH/oligo_files/oligos_$i.bed | cut -f 1 -d ' ')
    python3 ../util/compute_coverage_cv.py $DATA_PATH/aligned/exp_aligned_$i.ref $num_oligos $subsampling_coverage > $DATA_PATH/stats/exp_aligned_cvg_$subsampling_coverage.$i.cv.txt
done

# create .reads file for exp_1 (with unaligned reads)
sed -n '2~4p' $DATA_PATH/fastq/exp_1.fastq > $DATA_PATH/fastq/exp_1.reads
