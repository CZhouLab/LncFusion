# LncFusion

## Introduction

LncFusion is a computational pipeline designed to identify both lncRNA-involved fusions (lnc-fusions) and mRNA–mRNA fusions from RNA-seq data. To achieve this, we build a comprehensive index of lncRNA and mRNA by merging references from LncBook (v2) and GENCODE (v42). We then employ three fusion callers—STAR-Fusion, Arriba, and STAR-SEQR—to detect candidate fusions. Putative fusions are retained if they are called by at least two tools or meet an FFPM (fusion fragments per million) threshold ≥ 0.1 for STAR-Fusion-only calls. Mitochondrial, immunoglobulin, or highly duplicated genes are removed to reduce false positives. Finally, fusions are classified into lnc-fusions (lncRNA–lncRNA or lncRNA–mRNA) and mRNA-fusions (mRNA–mRNA)

![workflow]([Fig1.pdf](https://github.com/user-attachments/files/18494708/Fig1.pdf))

**Please cite our paper at *Non-Coding RNA* journal (https://doi.org/10.3390/ncrna8050070), if you find Flnc useful for your research. 

Version: 1.0.0

Last Modified: 07/27/2022

Authors: Zixiu Li (zixiu.li@umassmed.edu), Chan Zhou (chan.zhou@umassmed.edu)

Maintainer: Zixiu Li
