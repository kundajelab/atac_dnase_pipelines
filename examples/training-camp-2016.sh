#!/bin/bash

SUFFIX=Ct-1h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Ct-1h_S5_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Ct-1h_S5_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=Ct-3h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Ct-3h_S12_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Ct-3h_S12_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=Cu-1h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Cu-1h_S4_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Cu-1h_S4_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=Cu-3h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Cu-3h_S11_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Cu-3h_S11_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=Cz-1h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Cz-1h_S2_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Cz-1h_S2_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=Cz-3h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Cz-3h_S9_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Cz-3h_S9_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=DMSO-1h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/DMSO-1h_S6_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/DMSO-1h_S6_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=DMSO-3h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/DMSO-3h_S13_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/DMSO-3h_S13_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=Kz-1h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Kz-1h_S1_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Kz-1h_S1_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=Kz-3h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Kz-3h_S8_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Kz-3h_S8_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=Mz-1h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Mz-1h_S3_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Mz-1h_S3_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=Mz-3h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Mz-3h_S10_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/Mz-3h_S10_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=WT-1h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/WT-1h_S7_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/WT-1h_S7_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3

SUFFIX=WT-3h; WORK=/srv/scratch/shared/nandi/projects/training-camp-2016/run/$SUFFIX; mkdir -p $WORK; cd $WORK
FASTQ1_1=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/WT-3h_S14_L001_R1_001.fastq.gz
FASTQ1_2=/srv/scratch/shared/nandi/projects/training-camp-2016/data/fastqs/WT-3h_S14_L001_R2_001.fastq.gz
bds_scr $SUFFIX /users/leepc12/code/bds_atac/atac.bds -pe -title $SUFFIX -nth 5 -species saccer3 -fastq1_1 $FASTQ1_1 -fastq1_2 $FASTQ1_2 -url_base http://mitra.stanford.edu/kundaje/leepc12/training-camp-2016/$SUFFIX/out
sleep 3
