#!/bin/bash

ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds

BOWTIE_IDX=/srv/scratch/leepc12/hg19/bowtie2/ENCODEHg19_male
READ1=$1
READ2=$2
NUMTHREADS=4
GENOMESIZE=hs
CHROMSIZE=/mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes
V_INDEX=$3

bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod 'bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013' -shcmd 'export _JAVA_OPTIONS="-Xms256M -Xmx512M -XX:ParallelGCThreads=1"; export MAX_JAVA_MEM="4G"; export MALLOC_ARENA_MAX=4'

