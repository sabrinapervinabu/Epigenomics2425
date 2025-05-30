# ChIP-seq Analysis Pipeline
This repository contains a reproducible shell-based pipeline for analyzing ChIP-seq data, including quality filtering, peak calling, blacklist filtering, and comparison to ENCODE reference datasets.

## Overview
The pipeline performs the following steps:
1. **Download raw BAM files** (replicates and control).
2. **Filter uniquely mapped reads** using `samtools`.
3. **Generate alignment statistics**.
4. **Call peaks** with `MACS2` for individual replicates and merged data.
5. **Assess reproducibility** via peak intersections.
6. **Remove ENCODE blacklist regions** from peak sets.
7. **Compare filtered peaks** with ENCODE ChIP-seq reference peaks using the Jaccard index.

## Directory Structure
- `data/` – Input BAM files and blacklist/reference BED files
- `qc/` – Quality control and statistics
- `peaks/` – Output peak files from MACS2 and BEDTools
- `blacklist/` – Optional: Store blacklist files separately

## Run the pipeline 
```bash
chmod +x chipseq_pipeline.sh
./chipseq_pipeline.sh
```

