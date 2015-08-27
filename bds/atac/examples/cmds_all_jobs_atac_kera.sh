#!/bin.bash


# HEADER (hg19)
ATAC_BDS=/users/leepc12/code/pipelines/bds/atac/atac.bds
MOD_DEF="bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2; texlive/2013"
BOWTIE_IDX=/srv/scratch/leepc12/hg19/bowtie2/ENCODEHg19_male
NUMTHREADS=4
GENOMESIZE=hs
CHROMSIZE=/mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes
V_INDEX=/srv/scratch/csfoo/projects/atac_eugene/parsed_hg19_RefSeq.merged.ANS.bed

DATA_DIR=/srv/scratch/leepc12/data/atac_kera
RUN_DIR=/srv/scratch/leepc12/run/atac_kera



#BODY

SUBDIR=Froz_S11
READ1=${DATA_DIR}/Froz-Kera-Day-0-a-10min_S11_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-0-a-10min_S11_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S5
READ1=${DATA_DIR}/Froz-Kera-Day-0-a-1min_S5_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-0-a-1min_S5_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



SUBDIR=Froz_S12
READ1=${DATA_DIR}/Froz-Kera-Day-0-b-10min_S12_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-0-b-10min_S12_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S6
READ1=${DATA_DIR}/Froz-Kera-Day-0-b-1min_S6_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-0-b-1min_S6_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S13
READ1=${DATA_DIR}/Froz-Kera-Day-3-a-10min_S13_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-3-a-10min_S13_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S7
READ1=${DATA_DIR}/Froz-Kera-Day-3-a-1min_S7_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-3-a-1min_S7_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S14
READ1=${DATA_DIR}/Froz-Kera-Day-3-b-10min_S14_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-3-b-10min_S14_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S8
READ1=${DATA_DIR}/Froz-Kera-Day-3-b-1min_S8_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-3-b-1min_S8_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S15
READ1=${DATA_DIR}/Froz-Kera-Day-6-a-10min_S15_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-6-a-10min_S15_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S9
READ1=${DATA_DIR}/Froz-Kera-Day-6-a-1min_S9_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-6-a-1min_S9_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S16
READ1=${DATA_DIR}/Froz-Kera-Day-6-b-10min_S16_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-6-b-10min_S16_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=Froz_S10
READ1=${DATA_DIR}/Froz-Kera-Day-6-b-1min_S10_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Froz-Kera-Day-6-b-1min_S10_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"




SUBDIR=S9
READ1=${DATA_DIR}/Kera-Day-0-a-10min_S9_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-0-a-10min_S9_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S5
READ1=${DATA_DIR}/Kera-Day-0-a-1min_S5_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-0-a-1min_S5_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S10
READ1=${DATA_DIR}/Kera-Day-0-b-10min_S10_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-0-b-10min_S10_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S6
READ1=${DATA_DIR}/Kera-Day-0-b-1min_S6_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-0-b-1min_S6_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S11
READ1=${DATA_DIR}/Kera-Day-3-a-10min_S11_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-3-a-10min_S11_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S7
READ1=${DATA_DIR}/Kera-Day-3-a-1min_S7_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-3-a-1min_S7_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S12
READ1=${DATA_DIR}/Kera-Day-3-b-10min_S12_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-3-b-10min_S12_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S8
READ1=${DATA_DIR}/Kera-Day-3-b-1min_S8_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-3-b-1min_S8_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"




SUBDIR=S15
READ1=${DATA_DIR}/Kera-Day-6-a-10min_S15_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-6-a-10min_S15_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S13
READ1=${DATA_DIR}/Kera-Day-6-a-1min_S13_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-6-a-1min_S13_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S16
READ1=${DATA_DIR}/Kera-Day-6-b-10min_S16_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-6-b-10min_S16_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"


SUBDIR=S14
READ1=${DATA_DIR}/Kera-Day-6-b-1min_S14_L001_R1_001.fastq.gz
READ2=${DATA_DIR}/Kera-Day-6-b-1min_S14_L001_R2_001.fastq.gz
mkdir -p ${RUN_DIR}/$SUBDIR; cd ${RUN_DIR}/$SUBDIR; bds -s sge ${ATAC_BDS} ${BOWTIE_IDX} $READ1 $READ2 $NUMTHREADS $GENOMESIZE $CHROMSIZE ${V_INDEX} out -mod "${MOD_DEF}"



#REPORTING PDFs on MITRA
/users/leepc12/code/pipelines/bds/_tools/recursive_gather_ln.sh /srv/scratch/leepc12/run/atac_kera /srv/www/kundaje/leepc12/atac_kera_pdf report.pdf




