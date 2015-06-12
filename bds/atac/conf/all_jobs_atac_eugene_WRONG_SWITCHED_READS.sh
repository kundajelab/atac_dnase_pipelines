#!/bin/bash

echo DONT USE THIS, SOME READS ARE SWITCHED!

exit 1


#HUMAN

SUBDIR="3_Shallow/Human_CDC2_CGTACTAG_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC2_CGTACTAG_fastq/L6_CGTACTAG_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC2_CGTACTAG_fastq/L6_CGTACTAG_L006_R2_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_hg19_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="3_Shallow/Human_CDC1_AGGCAGAA_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC1_AGGCAGAA_fastq/L6_AGGCAGAA_L006_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC1_AGGCAGAA_fastq/L6_AGGCAGAA_L006_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_hg19_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="3_Shallow/Human_sirpa_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_sirpa_fastq/L6_CTCTCTAC_L006_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_sirpa_fastq/L6_CTCTCTAC_L006_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_hg19_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"


SUBDIR="3_Shallow/Human_CDC1_TAAGGCGA_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC1_TAAGGCGA_fastq/L6_TAAGGCGA_L006_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Human_CDC1_TAAGGCGA_fastq/L6_TAAGGCGA_L006_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_hg19_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"


#MOUSE

SUBDIR="3_Shallow/Brain_CAP_CD13neg_fastq_mouse_mm9"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Brain_CAP_CD13neg_fastq_mouse_mm9/L2_GCTACGCT_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Brain_CAP_CD13neg_fastq_mouse_mm9/L2_GCTACGCT_L002_R2_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="3_Shallow/SP_2_fastq_mouse"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/SP_2_fastq_mouse/L6_GGACTCCT_L006_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/SP_2_fastq_mouse/L6_GGACTCCT_L006_R2_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"


SUBDIR="3_Shallow/DP_2_fastq_mouse"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/DP_2_fastq_mouse/L6_TAGGCATG_L006_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/DP_2_fastq_mouse/L6_TAGGCATG_L006_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"


SUBDIR="3_Shallow/Brain_CAP_CD13pos_fastq_mouse_mm9"
READ1="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Brain_CAP_CD13pos_fastq_mouse_mm9/L2_CAGAGAGG_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/3_Shallow/Brain_CAP_CD13pos_fastq_mouse_mm9/L2_CAGAGAGG_L002_R2_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="1_Deep_Seq/DP_fastq_rep_mouse"
READ1="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/DP_fastq_rep_mouse/L6_TGGATCTG_L006_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/DP_fastq_rep_mouse/L6_TGGATCTG_L006_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

#NEW MOUSE
SUBDIR="1_Deep_Seq/PLN_LEC_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/PLN_LEC_fastq/L2_TTGACCCT_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/PLN_LEC_fastq/L2_TTGACCCT_L002_R2_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
#rm -rf /srv/scratch/leepc12/run/atac_eugene/$SUBDIR
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="1_Deep_Seq/PLN_FRC_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/PLN_FRC_fastq/L2_CCACTCCT_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/PLN_FRC_fastq/L2_CCACTCCT_L002_R2_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
#rm -rf /srv/scratch/leepc12/run/atac_eugene/$SUBDIR
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="2_Shallow_Seq_High_Priority/Astrocyte_142_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/Astrocyte_142_fastq/L2_GTGTGGTG_L002_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/Astrocyte_142_fastq/L2_GTGTGGTG_L002_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
#rm -rf /srv/scratch/leepc12/run/atac_eugene/$SUBDIR
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="2_Shallow_Seq_High_Priority/Microglia_140_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/Microglia_140_fastq/L2_AGGTTGGG_L002_R1_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/Microglia_140_fastq/L2_AGGTTGGG_L002_R2_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
#rm -rf /srv/scratch/leepc12/run/atac_eugene/$SUBDIR
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="2_Shallow_Seq_High_Priority/PLN_HEV_3-11_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/PLN_HEV_3-11_fastq/L2_TGCTGGGT_L002_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/PLN_HEV_3-11_fastq/L2_TGCTGGGT_L002_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
#rm -rf /srv/scratch/leepc12/run/atac_eugene/$SUBDIR
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="2_Shallow_Seq_High_Priority/PLN_CAP_3-11_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/PLN_CAP_3-11_fastq/L2_GAGGGGTT_L002_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/PLN_CAP_3-11_fastq/L2_GAGGGGTT_L002_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
#rm -rf /srv/scratch/leepc12/run/atac_eugene/$SUBDIR
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"

SUBDIR="2_Shallow_Seq_High_Priority/P_Neutro_CD11b+_fastq"
READ1="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/P_Neutro_CD11b+_fastq/L2_ACCACTGT_L002_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/2_Shallow_Seq_High_Priority/P_Neutro_CD11b+_fastq/L2_ACCACTGT_L002_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
#rm -rf /srv/scratch/leepc12/run/atac_eugene/$SUBDIR
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"


#FAILED

SUBDIR="1_Deep_Seq/SP_fastq_rep_mouse"
READ1="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/SP_fastq_rep_mouse/L6_ACCACTGT_L006_R2_001.fastq.gz"
READ2="/srv/scratch/leepc12/data/atac_eugene/1_Deep_Seq/SP_fastq_rep_mouse/L6_ACCACTGT_L006_R1_001.fastq.gz"
BDS=/users/leepc12/code/pipelines/bds/atac/conf/bds_mm9_leepc12_on_nandi_mitra.sh
mkdir -p /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; cd /srv/scratch/leepc12/run/atac_eugene/$SUBDIR; $BDS "$READ1" "$READ2"




