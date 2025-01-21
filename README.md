# LncFusion

## Introduction

LncFusion is a computational pipeline designed to identify both lncRNA-involved fusions (lnc-fusions) and mRNA–mRNA fusions from RNA-seq data. To achieve this, we build a comprehensive index of lncRNA and mRNA by merging references from LncBook (v2) and GENCODE (v42). We then employ three fusion callers—STAR-Fusion, Arriba, and STAR-SEQR—to detect candidate fusions. Putative fusions are retained if they are called by at least two tools or meet an FFPM (fusion fragments per million) threshold ≥ 0.1 for STAR-Fusion-only calls. Mitochondrial, immunoglobulin, or highly duplicated genes are removed to reduce false positives. Finally, fusions are classified into lnc-fusions (lncRNA–lncRNA or lncRNA–mRNA) and mRNA-fusions (mRNA–mRNA)

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

### Step 1: Download Flnc software from GitHub
```bash
cd /home/username/
git clone https://github.com/CZhouLab/Flnc
cd Flnc
chmod 755 *.py
```

### Step 2: Download LIB folder from Zenodo
```bash
cd /home/username/Flnc
wget -c https://zenodo.org/record/7853855/files/LIB.zip?download=1
unzip zenodo.org/record/7853855/files/LIB.zip\?download\=1
rm -f zenodo.org/record/7853855/files/LIB.zip\?download\=1
```

## Running Flnc

The Flnc tool has two subcommands single and pair. The single subcommand can take three types of input files: single-end RNA-seq data in FASTQ format, and transcript data either in BED format or in FASTA format. The pair subcommand can take two ends of the paired-end RNA-seq data in FASTQ format as the input.

**Note:**

-	Flnc can accept the FASTQ files uncompressed (.fastq) or compressed by gzip (.fastq.gz) as input files. 

-	The input file should be in Linux format. If the file was created in DOS/Windows, it should be converted to Linux format (e.g. using dos2unix. See https://phoenixnap.com/kb/convert-dos-to-unix for detail).  

-	The reference gene annotation should be the GTF format file for hg38 assembly.

-	**Make sure that singularity software have read permission for the input files, and have both read and write permission for files in the output folder. We recommend depositing the input files and specifying output folder to be on the same hard drive (disk) with the Flnc software.**

**Usage** 
```bash
python2 Flnc.py {pair,single} -l LIBRARY -o OUTPUT_DIR -f {fastq,fasta,bed} {-1 FILE1 -2 FILE2 | -u FILE} [optional options]

When running Flnc with paired RNA-seq data, it is critical that the *_1 files and the *_2 files of replicates appear in separate comma-delimited lists, and that the order of the files in the two lists is the same.

Subcommands:		choose one of the subcommands {pair,single}                

Arguments:

	-f, --format	The format of the input file: fastq, or fasta or bed.
                  	If using the pair subcommand, the format must be "fastq".	    
                  	If using single subcommand, the format can be fastq, or fasta, or bed.
	-1 FILE1	This argument is mandatory if using the pair subcommand. 
			Full path of the mate 1 file of paired FASTQ files, paired with the mate 2 file specified with "-2 " option.
			The mate 1 of replicates can be input through comma delimitation, e.g., "<path>/Rep1_1.fastq,<path>/Rep2_1.fastq".
	-2 FILE2	This argument is mandatory if using the pair subcommand.
			Full path of the mate 2 file of paired FASTQ files, paired with the mate 1 file specified with "-1 " option. 
			The mate 2 of replicates can be input through comma delimitation, e.g., "<path>/Rep1_2.fastq,<path>/Rep2_2.fastq".
	-u FILE		This argument is mandatory if using the single subcommand.
			Full path of the single input file. 
			If "-f fastq", please input the full path of FASTQ file of single-end RNA-seq data. FASTQ files for replicates can be input through comma delimitation, For example, "<path>/Rep1.fastq,<path>/Rep2.fastq". 
			If "-f fasta", please input the full path of files with transcripts in FASTA format.
			If "-f bed", please input the full path of files with transcripts in BED format.
	-l --library	Full path of the LIB folder, which can be downloaded from https://zenodo.org/record/7853855/files/LIB.zip?download=1
	-o --output_dir	Please specify the name of the output folder. This must be specified as a full path. For example, "-o /home/username/Flnc_sample1_output".

Options:
	-g --gtf_file	Full path of the reference gene annotation file in GTF format. 
			Default: gencode.v29.annotation.gtf in the LIB folder.
	-m --model	Choose the abbreviation of one of the following models: 
			rf: random forest
			lr: logistic regression
			nb: naïve Bayes
			dt: decision tree
			knn: k-nearest neighbors
			rbfsvm: support vector machines with RBF kernel
			lsvm: support vector machines with linear kernel
			ensemble: the common result predicted by all models 
			Default: rf
	-s --strand	This option is required only if "-f fastq", otherwise this argument is not needed.
			Specify strand-specific information with the following three options: 
			first: corresponds to fr-firststrand of the –library-type option in the TopHat tool for stranded RNA-seq data
			second: corresponds to fr-secondstrand of the –library-type option in the TopHat tool for stranded RNA-seq data
			unstrand: specific for unstranded RNA-seq data
			Default: first 
	-h/--help 	Show help message and exit
	-v/--version	Print version
```
