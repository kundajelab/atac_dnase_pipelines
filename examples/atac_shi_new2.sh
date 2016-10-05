#!/bin/bash

TITLE=A628T;SPECIES=hg19;
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160205_MONK_0465_BC8C8DACXX/L3/160205_MONK_0465_BC8C8DACXX_L3_TAAGGCGA_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160205_MONK_0465_BC8C8DACXX/L3/160205_MONK_0465_BC8C8DACXX_L3_TAAGGCGA_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A629P;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160205_MONK_0465_BC8C8DACXX/L3/160205_MONK_0465_BC8C8DACXX_L3_CGTACTAG_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160205_MONK_0465_BC8C8DACXX/L3/160205_MONK_0465_BC8C8DACXX_L3_CGTACTAG_2_pf.fastq.gz
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A62I3;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L1/160212_BRISCOE_0278_AC7YJJACXX_L1_CGTACTAG_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L1/160212_BRISCOE_0278_AC7YJJACXX_L1_CGTACTAG_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A62A6;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L1/160212_BRISCOE_0278_AC7YJJACXX_L1_AGGCAGAA_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L1/160212_BRISCOE_0278_AC7YJJACXX_L1_AGGCAGAA_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A62IK;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L1/160212_BRISCOE_0278_AC7YJJACXX_L1_TCCTGAGC_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L1/160212_BRISCOE_0278_AC7YJJACXX_L1_TCCTGAGC_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A6296;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L3/160212_BRISCOE_0278_AC7YJJACXX_L3_AGGCAGAA_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L3/160212_BRISCOE_0278_AC7YJJACXX_L3_AGGCAGAA_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A629G;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L3/160212_BRISCOE_0278_AC7YJJACXX_L3_TCCTGAGC_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L3/160212_BRISCOE_0278_AC7YJJACXX_L3_TCCTGAGC_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A629A;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L3/160212_BRISCOE_0278_AC7YJJACXX_L3_GGACTCCT_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L3/160212_BRISCOE_0278_AC7YJJACXX_L3_GGACTCCT_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A629O;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L3/160212_BRISCOE_0278_AC7YJJACXX_L3_TAGGCATG_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160212_BRISCOE_0278_AC7YJJACXX/L3/160212_BRISCOE_0278_AC7YJJACXX_L3_TAGGCATG_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A629Y;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160217_PINKERTON_0397_BC83FPACXX/L1/160217_PINKERTON_0397_BC83FPACXX_L1_TAAGGCGA_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160217_PINKERTON_0397_BC83FPACXX/L1/160217_PINKERTON_0397_BC83FPACXX_L1_TAAGGCGA_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A629Z;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160217_PINKERTON_0397_BC83FPACXX/L1/160217_PINKERTON_0397_BC83FPACXX_L1_AGGCAGAA_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160217_PINKERTON_0397_BC83FPACXX/L1/160217_PINKERTON_0397_BC83FPACXX_L1_AGGCAGAA_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A629W;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160217_PINKERTON_0397_BC83FPACXX/L1/160217_PINKERTON_0397_BC83FPACXX_L1_GGACTCCT_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160217_PINKERTON_0397_BC83FPACXX/L1/160217_PINKERTON_0397_BC83FPACXX_L1_GGACTCCT_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

TITLE=A62IE;SPECIES=hg19
FASTQ1_1=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160217_PINKERTON_0397_BC83FPACXX/L1/160217_PINKERTON_0397_BC83FPACXX_L1_TAGGCATG_1_pf.fastq.gz
FASTQ1_2=/srv/gsfs0/projects/seq_center/Illumina/PublishedResults/2016/feb/160217_PINKERTON_0397_BC83FPACXX/L1/160217_PINKERTON_0397_BC83FPACXX_L1_TAGGCATG_2_pf.fastq.gz
WORK=/srv/gsfs0/scratch/leepc12/run/atac_shi_new2/$TITLE; mkdir -p $WORK; cd $WORK
bds_scr $TITLE /home/leepc12/atac_dnase_pipelines/atac.bds -species $SPECIES -nth 8 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2

