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
from collections import defaultdict

def get_paired_end_read_count(r1_fa, r2_fa):
    cmd_r1 = [
        "apptainer",
        "exec",
        "Lib/seqkit.sif",
        "seqkit",
        "stats",
        "-T", r1_fa,
    ]
    cmd_r2 = [ 
        "apptainer",
        "exec",
        "Lib/seqkit.sif",
        "seqkit",
        "stats",
        "-T", r2_fa,
    ]
    result_r1 = subprocess.run(cmd_r1, capture_output=True, text=True, check=True)
    result_r2 = subprocess.run(cmd_r2, capture_output=True, text=True, check=True)

    num_seqs_r1 = int(result_r1.stdout.strip().split("\n")[1].split("\t")[3])
    num_seqs_r2 = int(result_r2.stdout.strip().split("\n")[1].split("\t")[3])

    if num_seqs_r1 != num_seqs_r2:
        raise ValueError(f"Read counts do not match: R1={count_r1}, R2={count_r2}")
    else:
        return num_seqs_r1

def STARFUSION(r1_fa, r2_fa, outputDir, CPU):
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir, exist_ok=True)
    
    cmd = [
        "apptainer",
        "exec",
        "Lib/LncFusion_v1.sif",
        "STAR-Fusion",
        "--left_fq", r1_fa,
        "--right_fq", r2_fa,
        "--genome_lib_dir", "Lib/LncFusion_v1_CTAT/ctat_genome_lib_build_dir",
        "--output_dir", outputDir,
        "--CPU", str(CPU)
    ]

    # Print the command (helpful for debugging/logging)
    logging.info("[INFO] Running STAR-Fusion with command:")
    logging.info("       " + " ".join(cmd))

    # Execute the command
    subprocess.run(cmd, check=True)
    os.system("rm -f " + outputDir + "/Aligned.out.bam")
    os.system("rm -f " + outputDir + "/Chimeric.out.junction")
    os.system("rm -f " + outputDir + "/ReadsPerGene.out.tab")
    os.system("rm -f " + outputDir + "/SJ.out.tab")
    os.system("rm -rf " + outputDir + "/_starF_checkpoints")
    os.system("rm -rf " + outputDir + "/star-fusion.preliminary")
    os.system("rm -rf " + outputDir + "/_STARgenome")
    os.system("rm -rf " + outputDir + "/_STARpass1")
    os.system("rm -rf " + outputDir + "/_STARtmp")
    TimeNow = str(datetime.datetime.now())
    logging.info("[INFO] STAR-Fusion completed successfully at %s", TimeNow)


def ARRIBA(r1_fa, r2_fa, outputDir, CPU):
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir, exist_ok=True)

    current_directory = os.getcwd()
    os.chdir(outputDir)
    
    cmd = [
        "apptainer",
        "exec",
        current_directory + "/Lib/LncFusion_v1.sif",
        "run_arriba.sh",
        current_directory + "/Lib/LncFusion_v1_ARRIBA_Reference/STAR_index_final/",
        current_directory + "/Lib/LncFusion_v1_ARRIBA_Reference/ref_annot.gtf",
        current_directory + "/Lib/LncFusion_v1_ARRIBA_Reference/ref_genome.fa",
        current_directory + "/Lib/LncFusion_v1_ARRIBA_Reference/blacklist_hg38_GRCh38_v2.3.0.tsv.gz",
        current_directory + "/Lib/LncFusion_v1_ARRIBA_Reference/known_fusions_hg38_GRCh38_v2.3.0.tsv.gz",
        current_directory + "/Lib/LncFusion_v1_ARRIBA_Reference/protein_domains_hg38_GRCh38_v2.3.0.gff3",
        str(CPU),    # Convert integer to string
        r1_fa,
        r2_fa
    ]

     # Print the command (helpful for debugging/logging)
    logging.info("[INFO] Running ARRIBA with command:")
    logging.info("       " + " ".join(cmd))

    # Execute the command
    subprocess.run(cmd, check=True)
    os.system("rm -f Aligned.sortedByCoord.out.bam")
    os.system("rm -f Aligned.sortedByCoord.out.bam.bai")
    os.system("rm -f SJ.out.tab")
    TimeNow = str(datetime.datetime.now())
    logging.info("[INFO] ARRIBA completed successfully at %s", TimeNow)
    os.chdir(current_directory)

def STARSEQR(r1_fa, r2_fa, outputDir, CPU):
    if os.path.exists(outputDir + "_STAR-SEQR"):
        os.system("rm -rf " + outputDir + "_STAR-SEQR")

    cmd = [
        "singularity",
        "exec",
        "Lib/starseqr_latest.sif",
        "starseqr.py",
        "-1", r1_fa,
        "-2", r2_fa,
        "-m", str(1),
        "-p", outputDir,
        "-t", str(CPU),
        "-i", "Lib/LncFusion_v1_STARSEQR_Index",
        "-g", "Lib/lncRNALncBookv2_GENCODEv42.gtf",
        "-r", "Lib/hg38.fa",
        "-vv"
    ]

    # Print the command (helpful for debugging/logging)
    logging.info("[INFO] Running STAR-SEQR with command:")
    logging.info("       " + " ".join(cmd))

    # Execute the command
    subprocess.run(cmd, check=True)
    TimeNow = str(datetime.datetime.now())
    logging.info("[INFO] STAR-SEQR completed successfully at %s", TimeNow)

def BP_detection(STARFusion_output, ARRIBA_output, STARSEQR_output, outputDir):
    STARFUSION_Select = {}
    STARFUSION_NoSelect = {}
    STARSEQR = {}
    ARRIBA = {}
    BP_DICT = {}
    BP_STARFUSION_JunctionReadCount = {}
    BP_STARFUSION_SpanningFragCount = {}
    BP_ARRIBA_JunctionReadCount = {}
    BP_ARRIBA_SpanningFragCount = {}
    BP_STARSEQR_JunctionReadCount = {}
    BP_STARSEQR_SpanningFragCount = {}

    f1 = open(STARFusion_output, "r")
    f1.readline()
    for line in f1:
        line = line.strip()
        element = line.split("\t")

        Left_Gene = element[0].split("-")[0]
        Left_BP = element[5]
        Right_Gene = element[0].split("-")[2]
        Right_BP = element[7]
        FFPM = element[11]

        BP_DICT[Left_BP+"_"+Right_BP] = 1
        BP_STARFUSION_JunctionReadCount[Left_BP+"_"+Right_BP] = element[1]
        BP_STARFUSION_SpanningFragCount[Left_BP+"_"+Right_BP] = element[2]

        if float(FFPM) > 0.1:
            STARFUSION_Select[Left_BP+"_"+Right_BP] = Left_Gene + "::" + Right_Gene
        else:
            STARFUSION_NoSelect[Left_BP+"_"+Right_BP] = Left_Gene + "::" + Right_Gene
    f1.close()

    f2 = open(STARSEQR_output, "r")
    f2.readline()
    for line in f2:
        line = line.strip()
        element = line.split("\t")

        Left_BP = element[6]
        Left_BP_chr = Left_BP.split(":")[0]
        Left_BP_pos = str(int(Left_BP.split(":")[1])+1)
        Left_BP_strand = Left_BP.split(":")[2]

        Right_BP = element[7]
        Right_BP_chr = Right_BP.split(":")[0]
        Right_BP_pos = str(int(Right_BP.split(":")[1])+1)
        Right_BP_strand = Right_BP.split(":")[2]

        BP_new = Left_BP_chr + ":" + Left_BP_pos + ":" + Left_BP_strand + "_" + Right_BP_chr + ":" + Right_BP_pos  + ":" + Right_BP_strand

        Left_Gene = element[0].split("-")[0]
        Right_Gene = element[0].split("-")[2]

        STARSEQR[BP_new] = Left_Gene + "::" + Right_Gene
        BP_DICT[BP_new] = 1

        BP_STARSEQR_JunctionReadCount[BP_new] = str(int(element[2]) + int(element[3]))
        BP_STARSEQR_SpanningFragCount[BP_new] = element[1]
    f2.close()

    f3 = open(ARRIBA_output, "r")
    f3.readline()
    for line in f3:
        line = line.strip()
        element = line.split("\t")

        Left_BP_strand = element[2].split("/")[0]
        Right_BP_strand = element[3].split("/")[0]

        Left_BP_pos = element[4]
        Right_BP_pos = element[5]

        BP_new = Left_BP_pos + ":" + Left_BP_strand + "_" + Right_BP_pos + ":" + Right_BP_strand

        Left_Gene = element[0]
        Right_Gene = element[1]

        ARRIBA[BP_new] = Left_Gene + "::" + Right_Gene
        BP_DICT[BP_new] = 1
        
        BP_ARRIBA_JunctionReadCount[BP_new] = str(int(element[9]) + int(element[10]))
        BP_ARRIBA_SpanningFragCount[BP_new] = element[11]
    f3.close() 

    output_file = outputDir + "/breakpoints_detection.txt"
    fo = open(output_file, "w")
    fo.write("Breakpoint\tToolsCount\tstarfusion_FFPM_high\tstarseqr\tarriba\tstarfusion_gene\tstarfusion_JunctionReadCount\tstarfusion_SpanningFragCount\tstarseqr_gene\tstarseqr_JunctionReadCount\tstarseqr_SpanningFragCount\tarriba_gene\tarriba_JunctionReadCount\tarriba_SpanningFragCount\n")
    for bp in BP_DICT:
        count = 0
        starfusion = ""
        starseqr = ""
        arriba = ""
        if bp in STARFUSION_Select:
            starfusion = "STARFUSION_Y"
            starfusion_gene = STARFUSION_Select[bp]
            starfusion_JunctionReadCount = BP_STARFUSION_JunctionReadCount[bp]
            starfusion_SpanningFragCount = BP_STARFUSION_SpanningFragCount[bp]
            count = count + 1
        elif bp in STARFUSION_NoSelect:
            starfusion = "STARFUSION_N"
            starfusion_gene = STARFUSION_NoSelect[bp]
            starfusion_JunctionReadCount = BP_STARFUSION_JunctionReadCount[bp]
            starfusion_SpanningFragCount = BP_STARFUSION_SpanningFragCount[bp]
            count = count + 1
        else:
            starfusion = "NA"
            starfusion_gene = "NA"
            starfusion_JunctionReadCount = 0
            starfusion_SpanningFragCount = 0

        if bp in STARSEQR:
            starseqr = "STARSEQR"
            starseqr_gene = STARSEQR[bp]
            starseqr_JunctionReadCount = BP_STARSEQR_JunctionReadCount[bp]
            starseqr_SpanningFragCount = BP_STARSEQR_SpanningFragCount[bp]
            count = count + 1
        else:
            starseqr = "NA"
            starseqr_gene = "NA"
            starseqr_JunctionReadCount = 0
            starseqr_SpanningFragCount = 0

        if bp in ARRIBA:
            arriba = "ARRIBA"
            arriba_gene = ARRIBA[bp]
            arriba_JunctionReadCount = BP_ARRIBA_JunctionReadCount[bp]
            arriba_SpanningFragCount = BP_ARRIBA_SpanningFragCount[bp]
            count = count + 1
        else:
            arriba = "NA"
            arriba_gene = "NA"
            arriba_JunctionReadCount = 0
            arriba_SpanningFragCount = 0

        fo.write(bp+"\t"+str(count)+"\t"+starfusion+"\t"+starseqr+"\t"+arriba+"\t"+starfusion_gene+"\t"+str(starfusion_JunctionReadCount)+"\t"+str(starfusion_SpanningFragCount)+"\t"+starseqr_gene+"\t"+str(starseqr_JunctionReadCount)+"\t"+str(starseqr_SpanningFragCount)+"\t"+arriba_gene+"\t"+str(arriba_JunctionReadCount)+"\t"+str(arriba_SpanningFragCount)+"\n")
    fo.close()

def Consensus_filtering(BP_detection_output, outputDir):
    Gene_Dict = {}
    STARFUSION_Y_Gene_Dict = {}
    STARFUSION_Gene_Dict = {}
    STARSEQR_Gene_Dict = {}
    ARRIBA_Gene_Dict = {}
    Gene_2_Dict = {}

    fi = open(BP_detection_output, "r")
    fi.readline()
    for line in fi:
        line = line.strip()
        element = line.split("\t")

        BP = element[0]
        ToolsCount = element[1]
        STARFUSION_Y = element[2]
        STARFUSION_Gene = element[5]
        STARSEQR_Gene = element[8]
        ARRIBA_Gene = element[11]
       
        Tools_list = list()
        if (int(ToolsCount) >= 2) or (STARFUSION_Y == "STARFUSION_Y"):
            if STARFUSION_Gene != "NA":
                Tools_list.append("STARFUSION")
            if STARSEQR_Gene != "NA":
                Tools_list.append("STARSEQR")
            if ARRIBA_Gene != "NA":
                Tools_list.append("ARRIBA")

            Tools_list_output = ";".join(Tools_list)
            if STARFUSION_Gene != "NA":
                Gene_2_Dict[STARFUSION_Gene] = STARFUSION_Y + "\t" + ToolsCount + "\t" + Tools_list_output
            if STARSEQR_Gene != "NA":
                Gene_2_Dict[STARSEQR_Gene] = STARFUSION_Y + "\t" + ToolsCount + "\t" + Tools_list_output
            if ARRIBA_Gene != "NA":
                Gene_2_Dict[ARRIBA_Gene] = STARFUSION_Y + "\t" + ToolsCount + "\t" + Tools_list_output

        if STARFUSION_Y == "STARFUSION_Y":
            STARFUSION_Y_Gene_Dict[STARFUSION_Gene] = 1
        else:
            pass

        if STARFUSION_Gene != "NA":
            STARFUSION_Gene_Dict[STARFUSION_Gene] = 1
            Gene_Dict[STARFUSION_Gene] = 1
        else:
            pass

        if STARSEQR_Gene != "NA":
            STARSEQR_Gene_Dict[STARSEQR_Gene] = 1
            Gene_Dict[STARSEQR_Gene] = 1
        else:
            pass

        if ARRIBA_Gene != "NA":
            ARRIBA_Gene_Dict[ARRIBA_Gene] = 1
            Gene_Dict[ARRIBA_Gene] = 1
        else:
            pass

    fi.close()

    fo = open(outputDir + "/Consensus_filtering.txt", "w")
    fo.write("FusionGene\tstarfusion_FFPM_high\tToolsCount\tTools\n")
    
    for gene in Gene_2_Dict:
        fo.write(gene+"\t"+Gene_2_Dict[gene]+"\n")

    for gene in Gene_Dict:
        if gene not in Gene_2_Dict:
            if gene in STARFUSION_Y_Gene_Dict:
                starfusion_y = "STARFUSION_Y"
            else:
                starfusion_y = "NA"

            Method_List = list()
            Method_Count = 0
            if gene in STARFUSION_Gene_Dict:
                Method_Count = Method_Count + 1
                Method_List.append("STARFUSION")
            else:
                pass

            if gene in STARSEQR_Gene_Dict:
                Method_Count = Method_Count + 1
                Method_List.append("STARSEQR")
            else:
                pass

            if gene in ARRIBA_Gene_Dict:
                Method_Count = Method_Count + 1
                Method_List.append("ARRIBA")
            else:
                pass

            method = ";".join(Method_List)
            if (starfusion_y == 1) or (Method_Count >= 2):
                fo.write(gene+"\t"+str(starfusion_y)+"\t"+str(Method_Count)+"\t"+method+"\n")
        
        fo.close()

def Other_filtering(Consensus_filtering_output, outputDir):
    GENE_Dict = {}
    f1 = open("Lib/Immunoglobulin_Genes.txt", "r")
    for line in f1:
        line = line.strip()
        element = line.split("\t")
        GENE = element[0]
        GENE_Dict[GENE] = 1
    f1.close()

    INDEX_DICT = defaultdict(dict)
    INDEX = 1
    f2 = open("Lib/STAR_FUSION_paralog_clusters.dat", "r")
    for line in f2:
        line = line.strip()
        element = line.split("\t")
        for e in element:
            INDEX_DICT[INDEX][e] = 1
        INDEX = INDEX + 1
    f2.close()

    fo = open(outputDir + "/Other_filtering.txt", "w")
    fi = open(Consensus_filtering_output, "r")
    title = fi.readline().strip()
    fo.write(title+"\n")
    for line in fi:
        line = line.strip()
        element = line.split("\t")
        FusionGene = element[0]
        gene_left = FusionGene.split(":")[0]
        gene_right = FusionGene.split(":")[2]

        Immu_label = "N"
        if (gene_left in GENE_Dict) or (gene_right in GENE_Dict):
            Immu_label = "Y"
        else:
            pass

        SameGene_label = "N"
        ParalogGene_label = "N"

        if gene_left == gene_right:
            SameGene_label = "Y"
        else:
            for index in INDEX_DICT:
                if (gene_left in INDEX_DICT[index]) and (gene_right in INDEX_DICT[index]):
                    ParalogGene_label = "Y"

        if (Immu_label == "N") and (SameGene_label == "N") and (ParalogGene_label == "N"):       
            fo.write(line+"\n")
    fi.close()
    fo.close()

def extend_genomic_coordinates(head_start, head_end, tail_start, tail_end, strand, extension=10000):
    if strand == "+":
        if head_end < tail_start:
            gap = tail_start - head_end
            direction = "intra-chromosomal,same strand,tail downstream"
        else:
            gap = head_start - tail_end
            direction = "intra-chromosomal,same strand,tail upstream"

    else:  # strand == "-"
        if head_start > tail_end:
            gap = head_start - tail_end
            direction = "intra-chromosomal,same strand,tail downstream"
        else:
            gap = tail_start - head_end
            direction = "intra-chromosomal,same strand,tail upstream"

    return {
        "gap": str(gap),
        "direction": direction
    }

def Annotation(Other_filtering_output, BP_detection_output, TotalReadCount, outputDir):
    Gene_Location = {}
    Gene_Biotype = {}
    FusionGene_Dict = {}
    BP_Gene_Dict = {}

    f1 = open("Lib/lncRNALncBookv2_GENCODEv42.gene.bed", "r")
    for line in f1:
        line = line.strip()
        element = line.split("\t")
        Chr = element[0]
        Start = element[1]
        End = element[2]
        Gene = element[3]
        Biotype = element[4]
        Strand = element[5]
        Gene_Location[Gene] = Chr + ":" + Start + "-" + End + ":" + Strand
        Gene_Biotype[Gene] = Biotype
    f1.close()

    fo1 = open(outputDir + "/output/FusionGene_Annotation.txt", "w")
    fo1.write("FusionGene\tFusionType\tHeadGene_Biotype\tTailGene_Biotype\tHeadGene_Location\tTailGene_Location\tGenomic_Oigination\tDistance(HeadGene,TailGene)\n")

    f2 = open(Other_filtering_output, "r")
    f2.readline()
    for line in f2:
        line = line.strip()
        element = line.split("\t")
        head_gene = element[0].split(":")[0]
        tail_gene = element[0].split(":")[2]
        FusionGene_Dict[element[0]] = 1

        head_gene_location = Gene_Location[head_gene]
        tail_gene_location = Gene_Location[tail_gene]
        
        head_gene_biotype = Gene_Biotype[head_gene]
        tail_gene_biotype = Gene_Biotype[tail_gene]
        if (head_gene_biotype == "protein_coding") and (tail_gene_biotype == "protein_coding"):
            fusion_type = "mRNA-fusion"
        elif (head_gene_biotype == "lncRNA") or (tail_gene_biotype == "lncRNA"):
            fusion_type = "lncRNA-fusion"
        else:
            fusion_type = "others"

        head_gene_chr = head_gene_location.split(":")[0]
        head_gene_start = head_gene_location.split(":")[1].split("-")[0]
        head_gene_end = head_gene_location.split(":")[1].split("-")[1]
        head_gene_strand = head_gene_location.split(":")[2]

        tail_gene_chr = tail_gene_location.split(":")[0]
        tail_gene_start = tail_gene_location.split(":")[1].split("-")[0]
        tail_gene_end = tail_gene_location.split(":")[1].split("-")[1]
        tail_gene_strand = tail_gene_location.split(":")[2]

        if (head_gene_chr != tail_gene_chr):
            fo1.write(element[0]+"\t"+fusion_type+"\t"+head_gene_biotype+"\t"+tail_gene_biotype+"\t"+head_gene_location+"\t"+tail_gene_location+"\tinter-chromosomal\tInf\n")
        elif (head_gene_chr == tail_gene_chr) and (head_gene_strand != tail_gene_strand):
            fo1.write(element[0]+"\t"+fusion_type+"\t"+head_gene_biotype+"\t"+tail_gene_biotype+"\t"+head_gene_location+"\t"+tail_gene_location+"\tintra-chromosomal,opposite strand\tInf\n")
        elif (head_gene_chr == tail_gene_chr) and (head_gene_strand == tail_gene_strand):
            extended_coordinates = extend_genomic_coordinates(
                int(head_gene_start), int(head_gene_end),
                int(tail_gene_start), int(tail_gene_end),
                head_gene_strand
            )
            fo1.write(element[0]+"\t"+fusion_type+"\t"+head_gene_biotype+"\t"+tail_gene_biotype+"\t"+head_gene_location+"\t"+tail_gene_location+"\t"+extended_coordinates['direction']+"\t"+extended_coordinates['gap']+"\n")

    f2.close()
    fo1.close()

    fo2 = open(outputDir + "/output/Breakpoint_Annotation.txt", "w")
    fo2.write("BreakPoint\tFusionGene\tFusionDetectionTool\tDistance(HeadBreakPoint,TailBreakPoint)\tJunctionReadCount\tSpanningFragCount\tFFPM\n")

    f3 = open(BP_detection_output, "r")
    f3.readline()
    for line in f3:
        line = line.strip()
        element = line.split("\t")
        BP = element[0]
        STARFUSION  = element[2]
        STARSEQR = element[3]
        ARRIBA = element[4]
        
        STARFUSION_Gene = element[5]
        starfusion_JunctionReadCount = int(element[6])
        starfusion_SpanningFragCount = int(element[7])
        STARSEQR_Gene = element[8]
        starseqr_JunctionReadCount = int(element[9])
        starseqr_SpanningFragCount = int(element[10])
        ARRIBA_Gene = element[11]
        arriba_JunctionReadCount = int(element[12])
        arriba_SpanningFragCount = int(element[13])

        Head_BP_Infor = BP.split("_")[0].split(":")
        Head_BP_Chr = Head_BP_Infor[0]
        Head_BP_Pos = Head_BP_Infor[1]
        Head_BP_Strand = Head_BP_Infor[2]

        Tail_BP_Infor = BP.split("_")[1].split(":")
        Tail_BP_Chr = Tail_BP_Infor[0]
        Tail_BP_Pos = Tail_BP_Infor[1]
        Tail_BP_Strand = Tail_BP_Infor[2]

        BP_Gap = "NA"
        if Head_BP_Chr != Tail_BP_Chr:
            BP_Gap = "Inf" 
        elif Head_BP_Strand != Tail_BP_Strand:
            BP_Gap = "Inf"
        else:
            BP_Gap = abs(int(Head_BP_Pos)-int(Tail_BP_Pos))

        JunctionReadCount = max(starfusion_JunctionReadCount, starseqr_JunctionReadCount, arriba_JunctionReadCount)
        SpanningFragCount = max(starfusion_SpanningFragCount, starseqr_SpanningFragCount, arriba_SpanningFragCount)
        FFPM = float(1000000 * (JunctionReadCount + SpanningFragCount) * TotalReadCount)

        gene_set = set()
        for gene in [STARFUSION_Gene, STARSEQR_Gene, ARRIBA_Gene]:
            if gene in FusionGene_Dict:
                gene_set.add(gene)

        method_list = list()
        if STARFUSION != "NA":
            method_list.append("STAR-Fusion")
        if STARSEQR != "NA":
            method_list.append("STAR-SEQR")
        if ARRIBA != "NA":
            method_list.append("Arriba")

        if len(gene_set) > 0:
            Gene_name = ";".join(gene_set)
            Method_name = ";".join(method_list)
            fo2.write(BP+"\t"+Gene_name+"\t"+Method_name+"\t"+str(BP_Gap)+"\t"+str(JunctionReadCount)+"\t"+str(SpanningFragCount)+"\t"+str(FFPM)+"\n")
    f3.close()

    fo2.close()

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
        R1.fq (Full path of left.fq file)
    ''')
)
required_group.add_argument(
    '-2', '--right_fq',
    type=str,
    required=True,
    help=textwrap.dedent('''
        R2.fq (Full path of right.fq file)
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

os.makedirs(os.path.join(output_dir, "output"), exist_ok=True)

cpu = args.CPU if args.CPU else 20

LOG_FILENAME = output_dir + "/logfile.log"
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)
logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)
logging.info("=========================================")

TotalReadCount = get_paired_end_read_count(args.left_fq, args.right_fq)

STARFUSION(args.left_fq, args.right_fq, output_dir+"/STARFUSION", cpu)

ARRIBA(args.left_fq, args.right_fq, output_dir+"/ARRIBA", cpu)

STARSEQR(args.left_fq, args.right_fq, output_dir+"/STARSEQR", cpu)

BP_detection(output_dir + "/STARFUSION/star-fusion.fusion_predictions.tsv", output_dir + "/ARRIBA/fusions.tsv", output_dir + "/STARSEQR_STAR-SEQR/STARSEQR_STAR-SEQR_breakpoints.txt", output_dir)

Consensus_filtering(output_dir + "/breakpoints_detection.txt", output_dir)

Other_filtering(output_dir + "/Consensus_filtering.txt", output_dir)

Annotation(output_dir+"/Other_filtering.txt", output_dir+"/breakpoints_detection.txt", TotalReadCount, output_dir)

logging.info("=========================================")
logging.info("Process completed successfully!")

