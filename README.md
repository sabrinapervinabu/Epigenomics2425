# ChIP-seq Analysis Pipeline
This repository contains a reproducible workflow for analyzing ChIP-seq data, combining a shell-based pipeline and R scripts for downstream analysis. It includes quality filtering, peak calling, blacklist filtering, chromatin state visualization, and gene annotation summary based on GREAT output. 

### Overview
The analysis pipeline is composed of three main parts:
## 1. **Bash pipeline (`chipseq_pipeline.sh`)**  
   Automates download, preprocessing, peak calling, blacklist filtering, and reproducibility checks using `samtools`, `MACS2`, and `BEDTools`.
  The pipeline performs the following steps:
  1. **Download raw BAM files** (replicates and control).
  2. **Filter uniquely mapped reads** using `samtools`.
  3. **Generate alignment statistics**.
  4. **Call peaks** with `MACS2` for individual replicates and merged data.
  5. **Assess reproducibility** via peak intersections.
  6. **Remove ENCODE blacklist regions** from peak sets.
  7. **Compare filtered peaks** with ENCODE ChIP-seq reference peaks using the Jaccard index.

   #### Directory Structure
   - `data/` – Input BAM files and blacklist/reference BED files
   - `qc/` – Quality control and statistics
   - `peaks/` – Output peak files from MACS2 and BEDTools
   - `blacklist/` – Optional: Store blacklist files separately
    
   #### Run the pipeline from the terminal
   ```bash
   chmod +x chipseq_pipeline.sh
   ./chipseq_pipeline.sh
   ```

## 2. **R script for chromatin state visualization (`chromatin_states_plot.R`)**
   `K562_ChromHMM_15states.bed` and `chromatin_states_count` are included in the repository to allow users to test chromatin state classification and plotting functionality.
   
   ```bash
    bedtools intersect -a merged_summits.bed -b K562_ChromHMM_15states.bed -wa -wb > summits_annotated.bed 
    awk '{print $9}' summits_annotated.bed | sort | uniq -c > chromatin_states_count
   ```
   
   `chromatin_states_plot.R` takes as input a tab-delimited text file (chromatin_states_count) containing two columns:
   - Counts: the number of peaks associated with each chromatin state.
   - Chromatin_States: the chromatin state label (e.g., Tss, Enh1, Tx, etc.).

These states are grouped into biologically meaningful categories such as TSS, Enhancers, Transcription, etc.
The script generates pie and bar plots showing the distribution of peaks across functional chromatin state categories.
   #### Run the pipeline from RStudio:
   ```r
   source("chromatin_states_plot.R")
   ```
## 3. **R script for GREAT results processing (`great_gene_count.R`)**  
   `great_gene_count.R` is designed to assist in post-processing of GREAT analysis results, specifically to answer the following biological questions related to transcription factor (TF) binding:

  - How many genes are associated with TF binding at the promoter region?
  - How many genes are associated with TF binding at distal regulatory regions?
This follows the convention:
  - Promoter = ±1 kb from the transcription start site (TSS)
  - Distal region = up to 30 kb from the TSS, but >1 kb away
  
  Given the output file from GREAT (based on the “summit” file of peaks), this script:
  - Extracts and splits gene-distance associations in the column Species.assembly..hg38
  - Cleans and parses numeric values from entries like GeneA(+132) or GeneB(-5400)
  - Classifies gene targets as:
  - Promoter-bound if distance is ≤ 1000 bp
  - Distal-bound if distance is > 1000 bp
  - Counts how many genes fall into each category
  
  Ensure the file great.txt is in your working directory.
  #### Run the pipeline from the terminal
  ```bash
     Rscript great_gene_count.R
  ```
  Results are printed to the console. 
  
  'great.txt' is included to support testing of the great_gene_count.R script for promoter vs. distal region gene assignment
