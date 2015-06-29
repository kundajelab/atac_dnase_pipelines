#!/usr/bin/env bash
set -o pipefail
set -o errexit
set -o nounset

BOWTIE_IDX=$1
READS=$2
NUMTHREADS=$3
VINDEX=$4

MAPQ_THRESH=30

module load bowtie/2.2.4
module load samtools/1.2
module load picard-tools/1.129
module load bedtools/2.23.0
module load python_anaconda/2.2.0

trim_galore ${READS}

# TODO: generalize to .gz and .fq.gz etc
TRIMMED_FASTQ="$(basename $READS .fastq)_trimmed.fq"
RAW_BAM_FILE="${TRIMMED_FASTQ/.fq/.aligned.bam}"
ALIGNLOG="${TRIMMED_FASTQ/.fq/.align.log}"
OFPREFIX="${TRIMMED_FASTQ/.fq/}"

bowtie2 -x "$BOWTIE_IDX" --threads "$NUMTHREADS" -U <(zcat -f "$TRIMMED_FASTQ") -S /dev/stdout 2> "$ALIGNLOG" | \
    samtools view -o "$RAW_BAM_FILE" -b /dev/stdin

SORTED_BAM_FILE="${RAW_BAM_FILE/.bam/.sort.bam}"

# preseq needs sorted input
samtools sort -Ttmp.${RAW_BAM_FILE} -l0 -Obam "$RAW_BAM_FILE" -o "$SORTED_BAM_FILE"
samtools index "$SORTED_BAM_FILE"

# generate enrichment plot
python makeTSSplot.py "$SORTED_BAM_FILE" "$VINDEX" 2000 "${TRIMMED_FASTQ/.fq/}" -q $MAPQ_THRESH

# compute % chrM
CHRM_QC="${TRIMMED_FASTQ/.fq/.chrm}"
samtools idxstats "$SORTED_BAM_FILE" | \
    awk '{ if ($1 == "chrM") chrM = $3; totalreads += $3 } END { print chrM/totalreads "\t(" chrM "/" totalreads ")" }' > "${CHRM_QC}"


# generate PRESEQ plot
>&2 echo "PRESEQ Analysis..."
# this script must be in PATH
plotPRESEQ.sh "$SORTED_BAM_FILE"

# =============================
# Remove unmapped, mate unmapped
# not primary alignment, reads failing platform
# Remove low MAPQ reads
# ==================
FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
MAPQ_THRESH=30

samtools view -F 1796 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} | \
    samtools sort -Ttmp.${RAW_BAM_FILE} -l0 -Obam /dev/stdin -o "${FILT_BAM_FILE}"

# ========================
# Mark duplicates
# ======================
>&2 echo "Mark duplicates"
TMP_FILT_BAM_FILE="${FILT_BAM_PREFIX}.dupmark.bam"
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"

JVM_OPTS="-Xmx4G"
java $JVM_OPTS -jar ${PICARDROOT}/picard.jar MarkDuplicates \
   INPUT="${FILT_BAM_FILE}" OUTPUT="${TMP_FILT_BAM_FILE}" \
   METRICS_FILE="${DUP_FILE_QC}" VALIDATION_STRINGENCY=LENIENT \
   ASSUME_SORTED=true REMOVE_DUPLICATES=false

mv "${TMP_FILT_BAM_FILE}" "${FILT_BAM_FILE}"

# ============================
# Remove duplicates
# Index final position sorted BAM
# ============================
>&2 echo "Remove duplicates"
FINAL_BAM_PREFIX="${OFPREFIX}.filt.nodup.srt"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

samtools view -F 1796 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}

# Index Final BAM file
samtools index ${FINAL_BAM_FILE} 
mv ${FINAL_BAM_FILE}.bai ${FINAL_BAM_INDEX_FILE}

samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}

# =============================
# Compute library complexity
# =============================
# sort by position and strand
# Obtain unique count statistics

PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

# PBC File output
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

bedtools bamtobed -i ${FILT_BAM_FILE} | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | \
    grep -v 'chrM' | \
    sort | \
    uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}

rm ${FILT_BAM_FILE}

echo "Postprocessing done - sorted final BAM file is:"
echo "${FINAL_BAM_FILE}"
