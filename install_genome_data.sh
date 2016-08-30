#!/bin/bash
# Stop on error
set -e

## species file name, conda environment name

SPECIES_FILE_BASENAME=bds_atac_species.conf
CONDA_ENV=bds_atac

## build index or not

BUILD_BWT2_IDX=1
BUILD_BWA_IDX=0

## show help

if [ "$#" -lt 2 ]; then
  echo
  echo "This script installs data for genome [GENOME] on a directory [DATA_DIR]."
  echo "Genome data files will be installed on [DATA_DIR]/[GENOME]."
  echo "A species file [DATA_DIR]/${SPECIES_FILE_BASENAME} will be generated and added to default.conf."
  echo
  echo "Supported genomes: hg19, mm9, hg38 (BETA), mm10 (BETA)"
  echo
  echo "Usage: ./install_genome_data.sh [GENOME] [DATA_DIR]"
  echo "  Example: ./install_genome_data.sh hg19 $TMP/genome_data"
  echo
  exit 0
fi

GENOME=$1
DATA_DIR=$2
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SPECIES_FILE=${DATA_DIR}/${SPECIES_FILE_BASENAME}

echo 

## show warning

if [[ $GENOME == "hg19" || $GENOME == "mm9" ]]; then
  echo
elif [[ $GENOME == "hg38" ]]; then
  echo "Warning: hg38 is BETA (GRCh38.p3). There is no ATAQC data, uniq. map. tracks and blacklist (IDR peaks will not be filtered)."
  echo "Press any key to continue..."
  read -n1
elif [[ $GENOME == "mm10" ]]; then
  echo "Warning: mm10 is BETA (GRCm38.p4). There is no ATAQC data, uniq. map. tracks and blacklist (IDR peaks will not be filtered)."
  echo "Press any key to continue..."
  read -n1
else
  echo "Error: unsupported genome $GENOME"
  exit 1
fi

echo 

## data URLs

if [ $GENOME == "hg19" ]; then

  CHRSZ="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.chrom.sizes"
  REF_FA="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz"
  SEQ="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes"
  SEQ_CHR_ARR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y) #chr$i.fa.gz
  GENSZ="hs"
  UMAP="https://personal.broadinstitute.org/anshul/projects/umap/encodeHg19Male/globalmap_k20tok54.tgz"
  BLACKLIST="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

  # data for ATAQC
  TSS_ENRICH="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/hg19_RefSeq_stranded.bed.gz"
  DNASE="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
  PROM="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
  ENH="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
  REG2MAP="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
  ROADMAP_META="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/eid_to_mnemonic.txt"

elif [ $GENOME == "mm9" ]; then

  CHRSZ="http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/encodeDCC/referenceSequences/male.mm9.chrom.sizes"
  REF_FA="http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/encodeDCC/referenceSequences/male.mm9.fa.gz"
  SEQ="http://hgdownload.cse.ucsc.edu/goldenpath/mm9/chromosomes"
  SEQ_CHR_ARR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 M X Y)
  GENSZ="mm"
  BLACKLIST="https://personal.broadinstitute.org/anshul/projects/mouse/blacklist/mm9-blacklist.bed.gz"

  # data for ATAQC
  TSS_ENRICH="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/mm9.gencode.vM1.transcript.TSS.bed"
  DNASE="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/mm9_dhs_universal_ucsc_v1.bed.gz"
  PROM="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/mm9_dhs_prom_ucsc_v1.bed.gz"
  ENH="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/mm9_dhs_enh_ucsc_v1.bed.gz"
  REG2MAP="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/dnase_avgs_merged_named.fseq.vals.gz"
  ROADMAP_META="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/accession_to_name.txt"

elif [ $GENOME == "hg38" ]; then

  CHRSZ="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
  REF_FA="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_23/GRCh38.p3.genome.fa.gz"
  SEQ="http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes"
  SEQ_CHR_ARR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y)
  GENSZ="hs"
  UMAP=""
  BLACKLIST=""

elif [ $GENOME == "mm9" ]; then

  CHRSZ="http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/encodeDCC/referenceSequences/male.mm9.chrom.sizes"
  REF_FA="http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/encodeDCC/referenceSequences/male.mm9.fa.gz"
  SEQ="http://hgdownload.cse.ucsc.edu/goldenpath/mm9/chromosomes"
  SEQ_CHR_ARR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 M X Y)
  GENSZ="mm"
  UMAP="https://personal.broadinstitute.org/anshul/projects/umap/mm9/globalmap_k20tok54.tgz"
  BLACKLIST="https://personal.broadinstitute.org/anshul/projects/mouse/blacklist/mm9-blacklist.bed.gz"

elif [ $GENOME == "mm10" ]; then

  CHRSZ="http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
  REF_FA="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M10/GRCm38.p4.genome.fa.gz"
  SEQ="http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes"
  SEQ_CHR_ARR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 M X Y)
  GENSZ="mm"
  UMAP=""
  BLACKLIST=""

fi

## get prefix of reference genome
REF_FA_PREFIX=$(basename ${REF_FA} .gz)

## create directories
mkdir -p ${DATA_DIR}/$GENOME
mkdir -p ${DATA_DIR}/$GENOME/seq

## download files
echo "Downloading files..."
cd ${DATA_DIR}/$GENOME
if [[ $UMAP != "" ]]; then wget -N -c $UMAP; fi
wget -N -c $CHRSZ
wget -N -c ${REF_FA}
if [[ $BLACKLIST != "" ]]; then wget -N -c $BLACKLIST; fi
mkdir -p ataqc && cd ataqc
if [[ $TSS_ENRICH != "" ]]; then wget -N -c $TSS_ENRICH; fi
if [[ $DNASE != "" ]]; then wget -N -c $DNASE; fi
if [[ $PROM != "" ]]; then wget -N -c $PROM; fi
if [[ $ENH != "" ]]; then wget -N -c $ENH; fi
if [[ $REG2MAP != "" ]]; then wget -N -c $REG2MAP; fi
if [[ $ROADMAP_META != "" ]]; then wget -N -c $ROADMAP_META; fi
cd ${DATA_DIR}/$GENOME
mkdir -p seq && cd seq
for i in ${SEQ_CHR_ARR[@]}; do wget -N -c $SEQ/chr$i.fa.gz; done

## extract files
echo "Extracting/processing data files..."
cd ${DATA_DIR}/$GENOME
if [ ! -f ${REF_FA_PREFIX} ]; then
  gzip -d -f -c ${REF_FA_PREFIX}.gz > ${REF_FA_PREFIX}
fi
if [[ $UMAP != "" ]]; then tar zxvf $(basename $UMAP) --skip-old-files; fi
for i in ${SEQ_CHR_ARR[@]}; do
  if [ ! -f seq/chr$i.fa ]; then
    gzip -f -d -c seq/chr$i.fa.gz > seq/chr$i.fa
  fi
done

## build index

if [ ${BUILD_BWT2_IDX} == 1 ]; then
  echo "Building bowtie2 index..."
  mkdir -p ${DATA_DIR}/$GENOME/bowtie2_index
  cd ${DATA_DIR}/$GENOME/bowtie2_index
  rm -f ${REF_FA_PREFIX}
  ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
  if [ ! -f ${REF_FA_PREFIX}.rev.1.bt2 ]; then
    source activate ${CONDA_ENV}
    bowtie2-build ${REF_FA_PREFIX} ${REF_FA_PREFIX}
  fi
fi

if [ ${BUILD_BWA_IDX} == 1 ]; then
  echo "Building bwa index..."
  mkdir -p ${DATA_DIR}/$GENOME/bwa_index
  cd ${DATA_DIR}/$GENOME/bwa_index
  rm -f ${REF_FA_PREFIX}
  ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
  if [ ! -f ${REF_FA_PREFIX}.sa ]; then
    source activate ${CONDA_ENV}
    bwa index ${REF_FA_PREFIX}
  fi
fi

## create species file

echo "Creating species file... (${SPECIES_FILE})"
cd ${DATA_DIR} && touch ${SPECIES_FILE}

if [[ $(grep "\[$GENOME\]" ${SPECIES_FILE} | wc -l) < 1 ]]; then
  if [[ $UMAP != "" ]]; then UMAP_PATH="${DATA_DIR}/$GENOME/$(basename $UMAP .tgz)"; fi
  if [[ $BLACKLIST != "" ]]; then BLACKLIST_PATH="${DATA_DIR}/$GENOME/$(basename $BLACKLIST)"; fi
  if [[ $TSS_ENRICH != "" ]]; then TSS_ENRICH_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $TSS_ENRICH)"; fi
  if [[ $DNASE != "" ]]; then DNASE_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $DNASE)"; fi
  if [[ $PROM != "" ]]; then PROM_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $PROM)"; fi
  if [[ $ENH != "" ]]; then ENH_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $ENH)"; fi
  if [[ $REG2MAP != "" ]]; then REG2MAP_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $REG2MAP)"; fi
  if [[ $ROADMAP_META != "" ]]; then ROADMAP_META_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $ROADMAP_META)"; fi

  echo -e "[$GENOME] # installed by install_genome_data.sh" >> ${SPECIES_FILE}
  echo -e "chrsz\t= ${DATA_DIR}/$GENOME/$(basename $CHRSZ)" >> ${SPECIES_FILE}
  echo -e "seq\t= ${DATA_DIR}/$GENOME/seq" >> ${SPECIES_FILE}
  echo -e "gensz\t= $GENSZ" >> ${SPECIES_FILE}
  if [[ $UMAP != "" ]]; then echo -e "umap\t= ${UMAP_PATH}" >> ${SPECIES_FILE}; fi
  if [ ${BUILD_BWT2_IDX} == 1 ]; then
    echo -e "bwt2_idx\t= ${DATA_DIR}/$GENOME/bowtie2_index/${REF_FA_PREFIX}" >> ${SPECIES_FILE}
  fi
  if [ ${BUILD_BWA_IDX} == 1 ]; then
    echo -e "bwa_idx\t= ${DATA_DIR}/$GENOME/bwa_index/${REF_FA_PREFIX}" >> ${SPECIES_FILE}
  fi
  echo -e "ref_fa\t= ${DATA_DIR}/$GENOME/${REF_FA_PREFIX}" >> ${SPECIES_FILE}
  if [[ ${BLACKLIST_PATH} != "" ]]; then echo -e "blacklist\t= ${BLACKLIST_PATH}" >> ${SPECIES_FILE}; fi

  if [[ ${TSS_ENRICH_PATH} != "" ]]; then echo -e "# data for ATAQC" >> ${SPECIES_FILE}; fi
  if [[ ${TSS_ENRICH_PATH} != "" ]]; then echo -e "tss_enrich\t= ${TSS_ENRICH_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${DNASE_PATH} != "" ]]; then echo -e "dnase\t= ${DNASE_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${PROM_PATH} != "" ]]; then echo -e "prom\t= ${PROM_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${ENH_PATH} != "" ]]; then echo -e "enh\t= ${ENH_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${REG2MAP_PATH} != "" ]]; then echo -e "reg2map\t= ${REG2MAP_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${ROADMAP_META_PATH} != "" ]]; then echo -e "roadmap_meta\t= ${ROADMAP_META_PATH}" >> ${SPECIES_FILE}; fi
  echo "" >> ${SPECIES_FILE}
fi

# Add species file to ./default.env
sed -i "/DEF_SPECIES_FILE/c\species_file\t= ${SPECIES_FILE}\t# DEF_SPECIES_FILE: DO NOT MODIFY THIS COMMENT (install_genome_data.sh WILL NOT WORK)." "${SCRIPT_DIR}/default.env"

echo "=== Genome data ($GENOME) has been successfully installed. ==="
