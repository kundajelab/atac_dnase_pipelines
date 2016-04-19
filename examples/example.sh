find . -type f -name '*.fq.gz' -exec bash -c 'mv "$1" "${1/.fq.gz/.fastq.gz}"' -- {} \;


WORK_ROOT_OLD=/srv/scratch/leepc12/run/atac_old/atac_Hardison
mkdir -p $WORK/out/align/rep1
mkdir -p $WORK/out/align/rep2
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep1/*.fastq.gz $WORK/out/align/rep1
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep1/*.bam $WORK/out/align/rep1
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep1/*.bai $WORK/out/align/rep1
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep1/*.qc $WORK/out/align/rep1
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep1/*.align.log $WORK/out/align/rep1
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep1/*.txt $WORK/out/align/rep1
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep2/*.fastq.gz $WORK/out/align/rep2
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep2/*.bam $WORK/out/align/rep2
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep2/*.bai $WORK/out/align/rep2
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep2/*.qc $WORK/out/align/rep2
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep2/*.align.log $WORK/out/align/rep2
cp -pr $WORK_ROOT_OLD/$SUFFIX/out/rep2/*.txt $WORK/out/align/rep2





screen -RD HDS_01

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR428BSK; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF455RDD.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF583HWE.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out

screen -RD HDS_02

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR280ZDP; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF440YYI.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF698NQU.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out

screen -RD HDS_03

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR498DQA; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF943WTV.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF788TPW.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out

screen -RD HDS_04

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR914PYX; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF692NRE.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF757AOM.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out

screen -RD HDS_05

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR793RAV; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF535AXY.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF637MAT.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out

screen -RD HDS_06

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR257PGU; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF133CWQ.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF248WQE.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out

screen -RD HDS_07

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR064IHX; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF463PIR.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF360YBL.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out

screen -RD HDS_08

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR136XSY; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF765DWT.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF194NYS.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out

screen -RD HDS_09

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR229QKB; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF343AEN.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF199CEQ.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out

screen -RD HDS_10

WORK_ROOT=/srv/scratch/leepc12/run/atac_Hardison; SUFFIX=ENCSR453AWR; WORK=$WORK_ROOT/$SUFFIX 
mkdir -p $WORK; cd $WORK; mkdir -p out
bds $CODE/bds_atac/atac.bds -species mm9 -se -nth_macs2 1 -nth_spp 1 -par_lvl 1 \
-fastq1 $DATA/atac_Hardison_SE/ENCFF621CAV.fastq.gz \
-fastq2 $DATA/atac_Hardison_SE/ENCFF217JOA.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Hardison/$SUFFIX/out








screen -RD Peng_01

WORK_ROOT=/srv/scratch/leepc12/run/atac_Peng; SUFFIX=Adult_mouse_cerebellum; WORK=$WORK_ROOT/$SUFFIX
mkdir -p $WORK; cd $WORK
bds $CODE/bds_atac/atac.bds -se -nth_spp 1 -nth_macs2 1 -par_lvl 1 -species mm9  \
-fastq1 /srv/scratch/leepc12/data/atac_Peng/Adult_mouse_cerebellum/rep1/15160_CTCTCTAC-GTAAGGAG_rep1_pooled.fastq.gz \
-fastq2 /srv/scratch/leepc12/data/atac_Peng/Adult_mouse_cerebellum/rep2/15160_CTCTCTAC-GTAAGGAG_rep2_pooled.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Peng/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

screen -RD Peng_02

WORK_ROOT=/srv/scratch/leepc12/run/atac_Peng; SUFFIX=e14.5forebrain; WORK=$WORK_ROOT/$SUFFIX
mkdir -p $WORK; cd $WORK
bds $CODE/bds_atac/atac.bds -se -nth_spp 1 -nth_macs2 1 -par_lvl 1 -species mm9  \
-fastq1 /srv/scratch/leepc12/data/atac_Peng/e14.5forebrain/rep1/14770_GGACTCCT-AGAGTAGA_pooled.fastq.gz \
-fastq2 /srv/scratch/leepc12/data/atac_Peng/e14.5forebrain/rep2/14770_GGACTCCT_pooled.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Peng/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

screen -RD Peng_03

WORK_ROOT=/srv/scratch/leepc12/run/atac_Peng; SUFFIX=e14.5forebrain2; WORK=$WORK_ROOT/$SUFFIX
mkdir -p $WORK; cd $WORK
bds $CODE/bds_atac/atac.bds -se -nth_spp 1 -nth_macs2 1 -par_lvl 1 -species mm9  \
-fastq1 /srv/scratch/leepc12/data/atac_Peng/e14.5forebrain2/rep1/14787_CGTACTAG-TAGATCGC_pooled.fastq.gz \
-fastq2 /srv/scratch/leepc12/data/atac_Peng/e14.5forebrain2/rep2/14787_CGTACTAG_pooled.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Peng/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

screen -RD Peng_04

WORK_ROOT=/srv/scratch/leepc12/run/atac_Peng; SUFFIX=C2C12_24h_ATAC; WORK=$WORK_ROOT/$SUFFIX
mkdir -p $WORK; cd $WORK
bds $CODE/bds_atac/atac.bds -se -nth_spp 1 -nth_macs2 1 -par_lvl 1 -species mm9 \
-fastq /srv/scratch/leepc12/data/atac_Peng/C2C12_24h_ATAC/14851_CGAGGCTG-TATCCTCT_pooled.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Peng/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

screen -RD Peng_05

WORK_ROOT=/srv/scratch/leepc12/run/atac_Peng; SUFFIX=C2C12_Exponential; WORK=$WORK_ROOT/$SUFFIX
mkdir -p $WORK; cd $WORK
bds $CODE/bds_atac/atac.bds -se -nth_spp 1 -nth_macs2 1 -par_lvl 1 -species mm9 \
-fastq /srv/scratch/leepc12/data/atac_Peng/C2C12_Exponential/16007_GTAGAGGA_pooled.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Peng/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

screen -RD Peng_06

WORK_ROOT=/srv/scratch/leepc12/run/atac_Peng; SUFFIX=C2C12_Exponential2; WORK=$WORK_ROOT/$SUFFIX
mkdir -p $WORK; cd $WORK
bds $CODE/bds_atac/atac.bds -se -nth_spp 1 -nth_macs2 1 -par_lvl 1 -species mm9 \
-fastq /srv/scratch/leepc12/data/atac_Peng/C2C12_Exponential2/16008_TAAGGCGA_pooled.fastq.gz \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_Peng/$SUFFIX/out -env $CODE/bds_atac/conf/test.env





screen -RD atac_shi_01

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Colon-Sigmoid/M; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151218_LYNLEY_0531_AC8480ACXX_L2_TAAGGCGA_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151218_LYNLEY_0531_AC8480ACXX_L2_TAAGGCGA_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151218_LYNLEY_0531_AC8480ACXX_L2_CGTACTAG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151218_LYNLEY_0531_AC8480ACXX_L2_CGTACTAG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 4 -nth_spp 1 -nth_macs2 1 -par_lvl 5 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

2

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Colon-Transverse/M; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151218_LYNLEY_0531_AC8480ACXX_L2_GGACTCCT_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151218_LYNLEY_0531_AC8480ACXX_L2_GGACTCCT_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151218_LYNLEY_0531_AC8480ACXX_L2_TAGGCATG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151218_LYNLEY_0531_AC8480ACXX_L2_TAGGCATG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

3

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Colon-Transverse/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_TAAGGCGA_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_TAAGGCGA_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_CGTACTAG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_CGTACTAG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

4

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Colon-Sigmoid/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_AGGCAGAA_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_AGGCAGAA_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_TCCTGAGC_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_TCCTGAGC_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

5

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=spleen/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_GGACTCCT_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_GGACTCCT_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_TAGGCATG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151202_BRISCOE_0269_AC83PEACXX_L3_TAGGCATG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

6

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Breast-Mammary/M; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151203_MARPLE_0344_BC8C86ACXX_L3_CGTACTAG_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151203_MARPLE_0344_BC8C86ACXX_L3_CGTACTAG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

7

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=stomach/M; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151203_MARPLE_0344_BC8C86ACXX_L3_AGGCAGAA_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151203_MARPLE_0344_BC8C86ACXX_L3_AGGCAGAA_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151203_MARPLE_0344_BC8C86ACXX_L3_TCCTGAGC_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151203_MARPLE_0344_BC8C86ACXX_L3_TCCTGAGC_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

8

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Breast-Mammary/F1; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_TAAGGCGA_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_TAAGGCGA_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_CGTACTAG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_CGTACTAG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

9

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=stomach/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_AGGCAGAA_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_AGGCAGAA_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

############# durga

screen -RD atac_shi_1

10

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Adipose-Subcutaneous/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_GGACTCCT_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_GGACTCCT_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_TAGGCATG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L1_TAGGCATG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

11

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Artery-Tibial/M; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L4_GGACTCCT_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L4_GGACTCCT_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L4_TAGGCATG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151208_PINKERTON_0393_AC88RLACXX_L4_TAGGCATG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

12

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=spleen/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151210_MONK_0462_AC8C29ACXX_L4_TAAGGCGA_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151210_MONK_0462_AC8C29ACXX_L4_TAAGGCGA_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151210_MONK_0462_AC8C29ACXX_L4_CGTACTAG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151210_MONK_0462_AC8C29ACXX_L4_CGTACTAG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

13

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Artery-Tibial/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151210_MONK_0462_AC8C29ACXX_L4_GGACTCCT_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151210_MONK_0462_AC8C29ACXX_L4_GGACTCCT_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151210_MONK_0462_AC8C29ACXX_L4_TAGGCATG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151210_MONK_0462_AC8C29ACXX_L4_TAGGCATG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

14

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Adipose-Visceral/M; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L1_CGTACTAG_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L1_CGTACTAG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

15

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Esophagus-Mucosa/M; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L1_TCCTGAGC_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L1_TCCTGAGC_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

16

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Esophagus-Gast.-Junction/M; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L1_TAGGCATG_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L1_TAGGCATG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

17

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Adipose-Visceral/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L2_CGTACTAG_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L2_CGTACTAG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

18

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Esophagus-Mucosa/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L2_TCCTGAGC_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L2_TCCTGAGC_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

19

WORK_ROOT="/srv/scratch/leepc12/run/atac_shi"; SUFFIX=Esophagus-Gast.-Junction/F; WORK=$WORK_ROOT/$SUFFIX; mkdir -p $WORK; cd $WORK;
FASTQ_1_1=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L2_GGACTCCT_1_pf.fastq.gz
FASTQ_1_2=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L2_GGACTCCT_2_pf.fastq.gz
FASTQ_2_1=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L2_TAGGCATG_1_pf.fastq.gz
FASTQ_2_2=$DATA/atac_shi/151215_BRISCOE_0271_AC8C2AACXX_L2_TAGGCATG_2_pf.fastq.gz
bds $CODE/bds_atac/atac.bds -species hg19 -nth_bwt2 1 -nth_spp 1 -nth_macs2 1 -par_lvl 1 \
-fastq1_1 $FASTQ_1_1 -fastq1_2 $FASTQ_1_2 -fastq2_1 $FASTQ_2_1 -fastq2_2 $FASTQ_2_2 \
-url_base http://mitra.stanford.edu/kundaje/leepc12/atac_shi/$SUFFIX/out -env $CODE/bds_atac/conf/test.env

