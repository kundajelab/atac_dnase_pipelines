#!/bin/bash

ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
ATAC_PATH=/users/leepc12/code/pipelines/atac

BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
READ1=$1
READ2=$2
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/leepc12/mm9/TSS_mouse.bed

bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod 'bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2' -shcmd 'export _JAVA_OPTIONS="-Xms256M -Xmx512M -XX:ParallelGCThreads=1"; export MAX_JAVA_MEM="4G"; export MALLOC_ARENA_MAX=4' -addpath ${ATAC_PATH}

