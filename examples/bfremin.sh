#!/bin/bash

TITLE=Butyrate
FASTQ1_1=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/But-1-60K-CAGAGAGG_S3_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/But-1-60K-CAGAGAGG_S3_R2_001.fastq.gz
FASTQ2_1=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/But-2-60K-AAGAGGCA_S4_R1_001.fastq.gz
FASTQ2_2=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/But-2-60K-AAGAGGCA_S4_R2_001.fastq.gz
WORKDIR=/srv/scratch/shared/surya/leepc12/run/bfremin/ATAC-SEQ/$TITLE; mkdir -p $WORKDIR; cd $WORKDIR
bds_scr $TITLE /users/leepc12/code/bds_atac/atac.bds -nth 12 -species hg19 -title $TITLE -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -fastq2_1 $FASTQ2_1 -fastq2_2 $FASTQ2_2 
sleep 2

TITLE=ButSulf
FASTQ1_1=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/ButSulf-1-50K-AGGTTGGG_S7_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/ButSulf-1-50K-AGGTTGGG_S7_R2_001.fastq.gz
FASTQ2_1=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/ButSulf-2-60K-TTGACCCT_S8_R1_001.fastq.gz
FASTQ2_2=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/ButSulf-2-60K-TTGACCCT_S8_R2_001.fastq.gz
WORKDIR=/srv/scratch/shared/surya/leepc12/run/bfremin/ATAC-SEQ/$TITLE; mkdir -p $WORKDIR; cd $WORKDIR
bds_scr $TITLE /users/leepc12/code/bds_atac/atac.bds -nth 12 -species hg19 -title $TITLE -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -fastq2_1 $FASTQ2_1 -fastq2_2 $FASTQ2_2 
sleep 2

TITLE=Control
FASTQ1_1=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/Control-1-60K-CGTACTAG_S1_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/Control-1-60K-CGTACTAG_S1_R2_001.fastq.gz
FASTQ2_1=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/Control-2-60K-GGACTCCT_S2_R1_001.fastq.gz
FASTQ2_2=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/Control-2-60K-GGACTCCT_S2_R2_001.fastq.gz
WORKDIR=/srv/scratch/shared/surya/leepc12/run/bfremin/ATAC-SEQ/$TITLE; mkdir -p $WORKDIR; cd $WORKDIR
bds_scr $TITLE /users/leepc12/code/bds_atac/atac.bds -nth 12 -species hg19 -title $TITLE -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -fastq2_1 $FASTQ2_1 -fastq2_2 $FASTQ2_2 
sleep 2

TITLE=Sulfide
FASTQ1_1=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/Sulf-1-60K-CGAGGCTG_S5_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/Sulf-1-60K-CGAGGCTG_S5_R2_001.fastq.gz
FASTQ2_1=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/Sulf-2-60K-TGCTGGGT_S6_R1_001.fastq.gz
FASTQ2_2=/srv/scratch/shared/surya/leepc12/data/bfremin/ATAC-SEQ/Sulf-2-60K-TGCTGGGT_S6_R2_001.fastq.gz
WORKDIR=/srv/scratch/shared/surya/leepc12/run/bfremin/ATAC-SEQ/$TITLE; mkdir -p $WORKDIR; cd $WORKDIR
bds_scr $TITLE /users/leepc12/code/bds_atac/atac.bds -nth 12 -species hg19 -title $TITLE -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -fastq2_1 $FASTQ2_1 -fastq2_2 $FASTQ2_2 
sleep 2


