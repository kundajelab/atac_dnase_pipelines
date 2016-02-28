#!/bin/bash

#HUMAN

####################################################
SUBDIR="3_Shallow/Human_CDC2_CGTACTAG_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC2_CGTACTAG_fastq/L6_CGTACTAG_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC2_CGTACTAG_fastq/L6_CGTACTAG_L006_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/hg19/bowtie2/ENCODEHg19_male
NUMTHREADS=4
GENOMESIZE=hs
CHROMSIZE=/mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/parsed_hg19_RefSeq.merged.ANS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


####################################################
SUBDIR="3_Shallow/Human_CDC1_AGGCAGAA_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC1_AGGCAGAA_fastq/L6_AGGCAGAA_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC1_AGGCAGAA_fastq/L6_AGGCAGAA_L006_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/hg19/bowtie2/ENCODEHg19_male
NUMTHREADS=4
GENOMESIZE=hs
CHROMSIZE=/mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/parsed_hg19_RefSeq.merged.ANS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="3_Shallow/Human_sirpa_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_sirpa_fastq/L6_CTCTCTAC_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_sirpa_fastq/L6_CTCTCTAC_L006_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/hg19/bowtie2/ENCODEHg19_male
NUMTHREADS=4
GENOMESIZE=hs
CHROMSIZE=/mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/parsed_hg19_RefSeq.merged.ANS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="3_Shallow/Human_CDC1_TAAGGCGA_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC1_TAAGGCGA_fastq/L6_TAAGGCGA_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC1_TAAGGCGA_fastq/L6_TAAGGCGA_L006_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/hg19/bowtie2/ENCODEHg19_male
NUMTHREADS=4
GENOMESIZE=hs
CHROMSIZE=/mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/parsed_hg19_RefSeq.merged.ANS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"




#MOUSE

####################################################
SUBDIR="3_Shallow/Brain_CAP_CD13neg_fastq_mouse_mm9"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Brain_CAP_CD13neg_fastq_mouse_mm9/L2_GCTACGCT_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Brain_CAP_CD13neg_fastq_mouse_mm9/L2_GCTACGCT_L002_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


####################################################
SUBDIR="3_Shallow/SP_2_fastq_mouse"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/SP_2_fastq_mouse/L6_GGACTCCT_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/SP_2_fastq_mouse/L6_GGACTCCT_L006_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="3_Shallow/DP_2_fastq_mouse"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/DP_2_fastq_mouse/L6_TAGGCATG_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/DP_2_fastq_mouse/L6_TAGGCATG_L006_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="3_Shallow/Brain_CAP_CD13pos_fastq_mouse_mm9"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Brain_CAP_CD13pos_fastq_mouse_mm9/L2_CAGAGAGG_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Brain_CAP_CD13pos_fastq_mouse_mm9/L2_CAGAGAGG_L002_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="1_Deep_Seq/DP_fastq_rep_mouse"
READ1="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/DP_fastq_rep_mouse/L6_TGGATCTG_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/DP_fastq_rep_mouse/L6_TGGATCTG_L006_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="1_Deep_Seq/PLN_LEC_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/PLN_LEC_fastq/L2_TTGACCCT_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/PLN_LEC_fastq/L2_TTGACCCT_L002_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="1_Deep_Seq/PLN_FRC_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/PLN_FRC_fastq/L2_CCACTCCT_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/PLN_FRC_fastq/L2_CCACTCCT_L002_R2_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="2_Shallow_Seq_High_Priority/Astrocyte_142_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/Astrocyte_142_fastq/L2_GTGTGGTG_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/Astrocyte_142_fastq/L2_GTGTGGTG_L002_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="2_Shallow_Seq_High_Priority/Microglia_140_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/Microglia_140_fastq/L2_AGGTTGGG_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/Microglia_140_fastq/L2_AGGTTGGG_L002_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="2_Shallow_Seq_High_Priority/PLN_HEV_3-11_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/PLN_HEV_3-11_fastq/L2_TGCTGGGT_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/PLN_HEV_3-11_fastq/L2_TGCTGGGT_L002_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="2_Shallow_Seq_High_Priority/PLN_CAP_3-11_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/PLN_CAP_3-11_fastq/L2_GAGGGGTT_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/PLN_CAP_3-11_fastq/L2_GAGGGGTT_L002_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="2_Shallow_Seq_High_Priority/P_Neutro_CD11b+_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/P_Neutro_CD11b+_fastq/L2_ACCACTGT_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/P_Neutro_CD11b+_fastq/L2_ACCACTGT_L002_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=4
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



####################################################
SUBDIR="1_Deep_Seq/SP_fastq_rep_mouse"
READ1="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/SP_fastq_rep_mouse/L6_ACCACTGT_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/SP_fastq_rep_mouse/L6_ACCACTGT_L006_R2_001.fastq.gz"
####################################################
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/mm9/Bowtie2Index/genome
NUMTHREADS=12
GENOMESIZE=mm
CHROMSIZE=/srv/scratch/leepc12/mm9/mm9.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/mm9.gencode.vM1.transcript.TSS.bed
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; 
bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"




#REPORTING PDFs on MITRA
rm -rf /srv/www/kundaje/leepc12/atac_pdf/*
kdir -p /srv/www/kundaje/leepc12/atac_pdf
/users/leepc12/code/pipelines/bds/_tools/recursive_gather_ln.sh /srv/scratch/leepc12/run/atac_eugene /srv/www/kundaje/leepc12/atac_pdf report.pdf





