// paired end?
PAIRED_END=1

// general
NTHREADS=4 // bwa (align_bwa), Rscript (xcor)
TMP="/srv/gsfs0/scratch/leepc12"

// bwa alignment
BWA_INDEX_NAME="/srv/gsfs0/projects/kundaje/commonRepository/indexes/bwa_indexes/encodeHg19Male/v0.7.10/encodeHg19Male_bwa-0.7.10.fa"
BWA_PARAM="-q 5 -l 32 -k 2"

// remove dup
MARKDUP="/srv/gs1/software/picard-tools/1.92/MarkDuplicates.jar"

// post alignment filtering
MAPQ_THRESH=30

// spp
RSCRIPT="/srv/gs1/software/R/R-2.15.1/bin/Rscript"
NREADS=15000000
NREADS_PER_MILLION=15
NPEAK=300000
SPEAK=220

