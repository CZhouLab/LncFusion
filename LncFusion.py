#!/usr/bin/env python
import os
import re
import subprocess
import pandas as pd
import numpy as np
import sys
import itertools
from copy import deepcopy
import copy
import logging
import argparse
import textwrap
import datetime
import time

def STARFUSION(r1_fa, r2_fa, outputDir, CPU):
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir, exist_ok=True)
    
    cmd = [
        "STAR-Fusion",
        "--left_fq", r1_fa,
        "--right_fq", r2_fa,
        "--genome_lib_dir", "LncFusion_Software/Lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir",
        "--output_dir", outputDir,
        "--CPU", CPU
    ]

    # Print the command (helpful for debugging/logging)
    logging.info("[INFO] Running STAR-Fusion with command:")
    logging.info("       " + " ".join(cmd))

    # Execute the command
    subprocess.run(cmd, check=True)
    TimeNow = str(datetime.datetime.now())
    logging.info("[INFO] STAR-Fusion completed successfully at %s", TimeNow)


def ARRIBA(r1_fa, r2_fa, outputDir, CPU):
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir, exist_ok=True)
    
    cmd = [
        "run_arriba.sh",
        "LncFusion_Software/Lib/LncFusion_v1_ARRIBA_Reference/STAR_index_final/",
        "LncFusion_Software/Lib/LncFusion_v1_ARRIBA_Reference/ref_annot.gtf",
        "LncFusion_Software/Lib/LncFusion_v1_ARRIBA_Reference/ref_genome.fa",
        "LncFusion_Software/Lib/LncFusion_v1_ARRIBA_Reference/blacklist_hg38_GRCh38_v2.3.0.tsv.gz",
        "LncFusion_Software/Lib/LncFusion_v1_ARRIBA_Reference/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz",
        "LncFusion_Software/Lib/LncFusion_v1_ARRIBA_Reference/protein_domains_hg38_GRCh38_Loci.gff3",
        str(CPU),    # Convert integer to string
        r1_fa,
        r2_fa
    ]

     # Print the command (helpful for debugging/logging)
    logging.info("[INFO] Running ARRIBA with command:")
    logging.info("       " + " ".join(cmd))

    # Execute the command
    subprocess.run(cmd, check=True)
    TimeNow = str(datetime.datetime.now())
    logging.info("[INFO] ARRIBA completed successfully at %s", TimeNow)

def STARSEQR(r1_fa, r2_fa, outputDir, CPU):
starseqr.py -1 r1_fa -2 r2_fa -m 1 -p outputDir -t CPU -i LncFusion_Software/Lib/LncFusion_v1_STARSEQR_Index -g /home/zixiu.li-umw/LncFusion/input/lncRNALncBookv2_GENCODEv42_GTF/lncRNALncBookv2_GENCODEv42.gtf -r ${genome_fa} -vv



# Define shared parser
parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="LncFusion: A method to identify lncRNA-derived fusion transcripts from RNA-seq data."
)

# Group for required arguments
required_group = parser.add_argument_group('Required arguments')
required_group.add_argument(
    '-1', '--left_fq',
    type=str,
    required=True,
    help=textwrap.dedent('''
        R1.fq (left.fq file)
    ''')
)
required_group.add_argument(
    '-2', '--right_fq',
    type=str,
    required=True,
    help=textwrap.dedent('''
        R2.fq (right.fq file)
    ''')
)

    # Group for optional arguments
optional_group = parser.add_argument_group('Optional arguments')
optional_group.add_argument(
    '-o', '--output_dir',
    type=str,
    default='./LncFusion_output',
    help=textwrap.dedent('''
        Directory where the output files will be saved.
        Default: ./LncFusion_output
    ''')
)
optional_group.add_argument(
    '-c', '--CPU',
    type=int,
    default=20,
    help=textwrap.dedent('''
        number of threads.
        Default: 20
    ''')
)

# Parsing arguments
args = parser.parse_args()

# Creating output directory if it does not exist
output_dir = args.output_dir if args.output_dir else "./LncFusion_output"
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

cpu = args.CPU if args.CPU else 20

