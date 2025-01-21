# LncFusion

## Introduction

LncFusion is a computational pipeline designed to identify both lncRNA-involved fusions (lnc-fusions) and mRNA–mRNA fusions from RNA-seq data. To achieve this, we build a comprehensive index of lncRNA and mRNA by merging references from LncBook (v2) and GENCODE (v42). We then employ three fusion callers—STAR-Fusion, Arriba, and STAR-SEQR—to detect candidate fusions. Putative fusions are retained if they are called by at least two tools or meet an FFPM (fusion fragments per million) threshold ≥ 0.1 for STAR-Fusion-only calls. Mitochondrial, immunoglobulin, or highly duplicated genes are removed to reduce false positives. Finally, fusions are classified into lnc-fusions (lncRNA–lncRNA or lncRNA–mRNA) and mRNA-fusions (mRNA–mRNA)

![workflow](Fig1.png)

**Please cite our paper at *medRxiv* (https://medrxiv.org/cgi/content/short/2025.01.16.25320696v1), if you find LncFusion useful for your research. 

Version: 1.0.0

Last Modified: 01/20/2025

Authors: Zixiu Li (zixiu.li@umassmed.edu), Chan Zhou (chan.zhou@umassmed.edu)

Maintainer: Zixiu Li

## Prerequisites

To use LncFusino, you will need the following programs in your PATH:

•       apptainer (>=1.3.5)

•       python3 

•       OS: high performance computing cluster in Linux (suggested)

•       Reference genome: hg38


## Installation

### Step 1: Download LncFusion software from GitHub
```bash
cd /home/username/
git clone https://github.com/CZhouLab/LncFusion.git
cd LncFusion
chmod 755 *.py
```

### Step 2: Download LIB folder from Zenodo (wait to update)


## Running LncFusion

**Note:**

-	LncFusion can accept the FASTQ files uncompressed (.fastq) or compressed by gzip (.fastq.gz) as input files. 

-	The reference gene annotation should be the GTF format file for hg38 assembly.

**Usage** 
```bash
LncFusion.py [-h] -1 LEFT_FQ -2 RIGHT_FQ [-o OUTPUT_DIR] [-c CPU]

LncFusion: A method to identify lncRNA-derived fusion transcripts from RNA-seq data.

Required arguments:

	-1 --left_fq 	Path of the mate 1 file of paired FASTQ files, paired with the mate 2 file specified with "-2 " option.
	-2 --right_fq 	Path of the mate 2 file of paired FASTQ files, paired with the mate 1 file specified with "-1 " option. 

Optional arguments:
	-o --output_dir	Path of the output folder.
			Default: ./LncFusion_output
	-c --CPU	number of threads.
			Default: 20
	-h/--help 	Show help message and exit
```
