#!/bin/bash
set -e  # Stop the script if any command fails

# Create directories to organize files
mkdir -p data results qc peaks blacklist

# Step 1: Download BAM files (replicates and control)
wget -O data/replicate1_unfiltered.bam https://www.encodeproject.org/files/ENCFF831RPZ/@@download/ENCFF831RPZ.bam
wget -O data/replicate2_unfiltered.bam https://www.encodeproject.org/files/ENCFF827KFN/@@download/ENCFF827KFN.bam
wget -O data/control_unfiltered.bam   https://www.encodeproject.org/files/ENCFF423WGF/@@download/ENCFF423WGF.bam

# Step 2: Filter for uniquely mapped reads (MAPQ >= 1)
samtools view -bq 1 data/replicate1_unfiltered.bam > data/replicate1_unique.bam
samtools view -bq 1 data/replicate2_unfiltered.bam > data/replicate2_unique.bam
samtools view -bq 1 data/control_unfiltered.bam   > data/control_unique.bam

# Step 3: Generate mapping statistics
samtools flagstat data/replicate1_unfiltered.bam > qc/rep1_unfiltered_flagstat.txt
samtools flagstat data/replicate1_unique.bam     > qc/rep1_unique_flagstat.txt
samtools flagstat data/replicate2_unfiltered.bam > qc/rep2_unfiltered_flagstat.txt
samtools flagstat data/replicate2_unique.bam     > qc/rep2_unique_flagstat.txt
samtools flagstat data/control_unfiltered.bam    > qc/control_unfiltered_flagstat.txt
samtools flagstat data/control_unique.bam        > qc/control_unique_flagstat.txt

# Step 4: Peak calling using MACS2
macs2 callpeak -t data/replicate1_unique.bam -c data/control_unfiltered.bam -g hs -n REP1 --outdir peaks
macs2 callpeak -t data/replicate2_unique.bam -c data/control_unfiltered.bam -g hs -n REP2 --outdir peaks

# Step 5: Count number of peaks
wc -l peaks/REP1_peaks.narrowPeak
wc -l peaks/REP2_peaks.narrowPeak

# Step 6: Peak intersection between replicates
bedtools intersect -a peaks/REP1_peaks.narrowPeak -b peaks/REP2_peaks.narrowPeak | wc -l
bedtools intersect -a peaks/REP2_peaks.narrowPeak -b peaks/REP1_peaks.narrowPeak | wc -l

# Step 7: Call peaks using both replicates
macs2 callpeak -t data/replicate1_unique.bam data/replicate2_unique.bam -c data/control_unfiltered.bam -g hs -n merged --outdir peaks
wc -l peaks/merged_peaks.narrowPeak

# Step 8: Compare merged peaks with individual replicates
bedtools intersect -a peaks/REP1_peaks.narrowPeak -b peaks/merged_peaks.narrowPeak -u | wc -l
bedtools intersect -a peaks/REP2_peaks.narrowPeak -b peaks/merged_peaks.narrowPeak -u | wc -l

# Step 9: Filter peaks by ENCODE blacklist
wget -O data/blacklist.bed.gz https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
gunzip -f data/blacklist.bed.gz

bedtools intersect -v -a peaks/merged_peaks.narrowPeak -b data/blacklist.bed > peaks/merged_filtered_peaks.narrowPeak
wc -l peaks/merged_filtered_peaks.narrowPeak

bedtools intersect -v -a peaks/REP1_peaks.narrowPeak -b data/blacklist.bed > peaks/REP1_filtered_peaks.narrowPeak
bedtools intersect -v -a peaks/REP2_peaks.narrowPeak -b data/blacklist.bed > peaks/REP2_filtered_peaks.narrowPeak

bedtools intersect -a peaks/REP1_peaks.narrowPeak -b peaks/REP2_peaks.narrowPeak -u > peaks/intersect_peaks.narrowPeak
bedtools intersect -v -a peaks/intersect_peaks.narrowPeak -b data/blacklist.bed > peaks/intersect_filtered_peaks.narrowPeak
wc -l peaks/intersect_filtered_peaks.narrowPeak

# Step 10: Jaccard index with ENCODE reference peaks
wget -O data/EncodeNarrowPeak.bed.gz https://www.encodeproject.org/files/ENCFF310CCS/@@download/ENCFF310CCS.bed.gz
gunzip -f data/EncodeNarrowPeak.bed.gz
sort -k1,1 -k2,2n -k3,3n data/EncodeNarrowPeak.bed > data/EncodeNarrowPeak_sorted.bed

bedtools jaccard -a peaks/intersect_filtered_peaks.narrowPeak -b data/EncodeNarrowPeak_sorted.bed > qc/jaccard_intersect.txt
bedtools jaccard -a peaks/merged_filtered_peaks.narrowPeak -b data/EncodeNarrowPeak_sorted.bed > qc/jaccard_merged.txt
