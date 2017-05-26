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
  echo "Supported genomes: hg19, mm9, hg38 and mm10"
  echo
  echo "Usage: ./install_genome_data.sh [GENOME] [DATA_DIR]"
  echo "  Example: ./install_genome_data.sh hg19 $TMP/genome_data"
  echo
  exit 0
fi

GENOME=$1
DATA_DIR=$(readlink -f $2)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SPECIES_FILE=${DATA_DIR}/${SPECIES_FILE_BASENAME}
echo 

## data URLs

if [ $GENOME == "hg19" ]; then

  REF_FA="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz"
  UMAP="http://mitra.stanford.edu/kundaje/genome_data/hg19/globalmap_k20tok54.tgz"
  BLACKLIST="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz"

  # data for ATAQC
  TSS_ENRICH="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/hg19_gencode_tss_unique.bed.gz"
  DNASE="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.bed.gz"
  PROM="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/reg2map_honeybadger2_dnase_prom_p2.bed.gz"
  ENH="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/reg2map_honeybadger2_dnase_enh_p2.bed.gz"
  REG2MAP="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/dnase_avgs_reg2map_p10_merged_named.pvals.gz"
  ROADMAP_META="http://mitra.stanford.edu/kundaje/genome_data/hg19/ataqc/eid_to_mnemonic.txt"

elif [ $GENOME == "mm9" ]; then

  REF_FA="http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit"
  UMAP="http://mitra.stanford.edu/kundaje/genome_data/mm9/globalmap_k20tok54.tgz"
  BLACKLIST="http://mitra.stanford.edu/kundaje/genome_data/mm9/mm9-blacklist.bed.gz"

  # data for ATAQC
  TSS_ENRICH="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/mm9_gencode_tss_unique.bed.gz"
  DNASE="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/mm9_univ_dhs_ucsc.from_mm10.bed.gz"
  PROM="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/tss_mm9_master.from_mm10.bed.gz"
  ENH="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/mm9_enh_dhs_ucsc.from_mm10.bed.gz"
  REG2MAP_BED="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/mm9_dhs_universal_ucsc_v1.bed.gz"
  REG2MAP="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/dnase_avgs_merged_named.fseq.vals.gz"
  ROADMAP_META="http://mitra.stanford.edu/kundaje/genome_data/mm9/ataqc/accession_to_name.txt"

elif [ $GENOME == "hg38" ]; then

  REF_FA="https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
  BLACKLIST="http://mitra.stanford.edu/kundaje/genome_data/hg38/hg38.blacklist.bed.gz"

  # data for ATAQC
  TSS_ENRICH="http://mitra.stanford.edu/kundaje/genome_data/hg38/ataqc/hg38_gencode_tss_unique.bed.gz"
  DNASE="http://mitra.stanford.edu/kundaje/genome_data/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz"
  PROM="http://mitra.stanford.edu/kundaje/genome_data/hg38/ataqc/reg2map_honeybadger2_dnase_prom_p2.hg19_to_hg38.bed.gz"
  ENH="http://mitra.stanford.edu/kundaje/genome_data/hg38/ataqc/reg2map_honeybadger2_dnase_enh_p2.hg19_to_hg38.bed.gz"
  REG2MAP_BED="http://mitra.stanford.edu/kundaje/genome_data/hg38/ataqc/hg38_celltype_compare_subsample.bed.gz"
  REG2MAP="http://mitra.stanford.edu/kundaje/genome_data/hg38/ataqc/hg38_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="http://mitra.stanford.edu/kundaje/genome_data/hg38/ataqc/hg38_dnase_avg_fseq_signal_metadata.txt"

elif [ $GENOME == "mm10" ]; then

  REF_FA="https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz"
  BLACKLIST="http://mitra.stanford.edu/kundaje/genome_data/mm10/mm10.blacklist.bed.gz"

  # data for ATAQC
  TSS_ENRICH="http://mitra.stanford.edu/kundaje/genome_data/mm10/ataqc/mm10_gencode_tss_unique.bed.gz"
  DNASE="http://mitra.stanford.edu/kundaje/genome_data/mm10/ataqc/mm10_univ_dhs_ucsc.bed.gz"
  PROM="http://mitra.stanford.edu/kundaje/genome_data/mm10/ataqc/tss_mm10_master.bed.gz"
  ENH="http://mitra.stanford.edu/kundaje/genome_data/mm10/ataqc/mm10_enh_dhs_ucsc.bed.gz"
  REG2MAP_BED="http://mitra.stanford.edu/kundaje/genome_data/mm10/ataqc/mm10_celltype_compare_subsample.bed.gz"
  REG2MAP="http://mitra.stanford.edu/kundaje/genome_data/mm10/ataqc/mm10_dnase_avg_fseq_signal_formatted.txt.gz"
  ROADMAP_META="http://mitra.stanford.edu/kundaje/genome_data/mm10/ataqc/mm10_dnase_avg_fseq_signal_metadata.txt"

elif [ $GENOME == "macam7" ]; then

  REF_FA="http://www.unmc.edu/rhesusgenechip/MacaM_Rhesus_Genome_v7.fasta.bz2"
  EXTRA_LINE="nonamecheck = true # for bedtools >= 2.24. this prevents name convention error in bedtools intersect"

elif [ $GENOME == "dm3" ]; then
  
  REF_FA="http://hgdownload-test.cse.ucsc.edu/goldenPath/dm3/bigZips/dm3.2bit"

elif [ $GENOME == "pantro5" ]; then

  REF_FA="http://hgdownload.cse.ucsc.edu/goldenPath/panTro5/bigZips/panTro5.fa.gz"
fi

## show warning

if [[ ${REF_FA} == "" ]]; then
  echo "Error: unsupported genome $GENOME"
  exit 1
fi

if [[ ${BLACKLIST} == "" ]]; then
  echo "Warning: blacklist is not provided for $GENOME. IDR peaks will not be filtered."
  echo "Press any key to continue..."
  read -n1
fi

if [[ ${UMAP} == "" ]] && [[ ${CONDA_ENV} == "aquas_chipseq" ]]; then
  echo "Warning: unique mappability tracks are not provided for $GENOME."
  echo "Press any key to continue..."
  read -n1
fi

if [[ ${TSS_ENRICH} == "" ]] && [[ ${CONDA_ENV} == "bds_atac" ]]; then
  echo "Warning: Data for ATAQC are not provided for $GENOME. No ATAQC will be available."
  echo "Press any key to continue..."
  read -n1
fi

echo 

## create directories
mkdir -p ${DATA_DIR}/$GENOME
mkdir -p ${DATA_DIR}/$GENOME/seq

## download files
echo "Downloading files..."
cd ${DATA_DIR}/$GENOME
if [[ $UMAP != "" ]]; then wget -N -c $UMAP; fi
wget -c -O $(basename ${REF_FA}) ${REF_FA}
if [[ $BLACKLIST != "" ]]; then wget -N -c $BLACKLIST; fi
mkdir -p ataqc && cd ataqc
if [[ $TSS_ENRICH != "" ]]; then wget -N -c $TSS_ENRICH; fi
if [[ $DNASE != "" ]]; then wget -N -c $DNASE; fi
if [[ $PROM != "" ]]; then wget -N -c $PROM; fi
if [[ $ENH != "" ]]; then wget -N -c $ENH; fi
if [[ $REG2MAP != "" ]]; then wget -N -c $REG2MAP; fi
if [[ $REG2MAP_BED != "" ]]; then wget -N -c $REG2MAP_BED; fi
if [[ $ROADMAP_META != "" ]]; then wget -N -c $ROADMAP_META; fi

## extract unique mappability tracks
cd ${DATA_DIR}/$GENOME
if [[ $UMAP != "" ]]; then tar zxvf $(basename $UMAP) --skip-old-files; fi

## extract fasta and get prefix of reference genome
echo "Extracting/processing data files..."
cd ${DATA_DIR}/$GENOME

source activate ${CONDA_ENV}

if [[ ${REF_FA} == *.gz ]]; then 
  REF_FA_PREFIX=$(basename ${REF_FA} .gz)
  gzip -d -f -c ${REF_FA_PREFIX}.gz > ${REF_FA_PREFIX}
elif [[ ${REF_FA} == *.bz2 ]]; then
  REF_FA_PREFIX=$(basename ${REF_FA} .bz2)
  bzip2 -d -f -c ${REF_FA_PREFIX}.bz2 > ${REF_FA_PREFIX}
elif [[ ${REF_FA} == *.2bit ]]; then
  REF_FA_PREFIX=$(basename ${REF_FA} .2bit).fa
  twoBitToFa $(basename ${REF_FA}) ${REF_FA_PREFIX}
else
  REF_FA_PREFIX=$(basename ${REF_FA})  
fi

## extract fasta per chromosome
cd ${DATA_DIR}/$GENOME
mkdir -p seq
cd seq
rm -f ${REF_FA_PREFIX}
ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
faidx -x ${REF_FA_PREFIX}
cp --remove-destination *.fai ../

## create chrom sizes file
CHRSZ=$GENOME.chrom.sizes
#cut -f1,2 ${REF_FA_PREFIX}.fai | grep chr > ../$CHRSZ
cut -f1,2 ${REF_FA_PREFIX}.fai > ../$CHRSZ

## determine gensz
cd ${DATA_DIR}/$GENOME
GENSZ=$(cat $CHRSZ | awk '{sum+=$2} END{print sum}')
if [[ $GENOME == hg* ]]; then GENSZ=hs; fi
if [[ $GENOME == mm* ]]; then GENSZ=mm; fi

## build index
if [ ${BUILD_BWT2_IDX} == 1 ]; then
  echo "Building bowtie2 index..."
  mkdir -p ${DATA_DIR}/$GENOME/bowtie2_index
  cd ${DATA_DIR}/$GENOME/bowtie2_index
  rm -f ${REF_FA_PREFIX}
  ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
  if [ ! -f ${REF_FA_PREFIX}.rev.1.bt2 ]; then
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
  if [[ $REG2MAP_BED != "" ]]; then REG2MAP_BED_PATH="${DATA_DIR}/$GENOME/ataqc/$(basename $REG2MAP_BED)"; fi
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
  if [[ ${SPECIES_BROWSER} != "" ]]; then echo -e "species_browser\t= ${SPECIES_BROWSER}" >> ${SPECIES_FILE}; fi

  if [[ ${TSS_ENRICH_PATH} != "" ]]; then echo -e "# data for ATAQC" >> ${SPECIES_FILE}; fi
  if [[ ${TSS_ENRICH_PATH} != "" ]]; then echo -e "tss_enrich\t= ${TSS_ENRICH_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${DNASE_PATH} != "" ]]; then echo -e "dnase\t= ${DNASE_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${PROM_PATH} != "" ]]; then echo -e "prom\t= ${PROM_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${ENH_PATH} != "" ]]; then echo -e "enh\t= ${ENH_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${REG2MAP_PATH} != "" ]]; then echo -e "reg2map\t= ${REG2MAP_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${REG2MAP_BED_PATH} != "" ]]; then echo -e "reg2map_bed\t= ${REG2MAP_BED_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${ROADMAP_META_PATH} != "" ]]; then echo -e "roadmap_meta\t= ${ROADMAP_META_PATH}" >> ${SPECIES_FILE}; fi
  if [[ ${EXTRA_LINE} != "" ]]; then echo -e "${EXTRA_LINE}" >> ${SPECIES_FILE}; fi
  echo "" >> ${SPECIES_FILE}
fi

# Add species file to ./default.env
sed -i "/DEF_SPECIES_FILE/c\species_file\t= ${SPECIES_FILE}\t# DEF_SPECIES_FILE: DO NOT MODIFY THIS COMMENT (install_genome_data.sh WILL NOT WORK)." "${SCRIPT_DIR}/default.env"

echo "=== Genome data ($GENOME) has been successfully installed. ==="
