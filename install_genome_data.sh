#!/bin/bash
# Stop on error
set -e

SPECIES_FILE_BASENAME=bds_atac_species.conf
CONDA_ENV=bds_atac

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

if [[ $GENOME == "hg19" || $GENOME == "mm9" ]]; then
  echo
elif [[ $GENOME == "hg38" ]]; then
  echo "Warning: hg38 is BETA (GRCh38.p3). There is no ATAQC data and blacklist (IDR peaks will not be filtered)."
  echo "Press any key to continue..."
  read -n1
elif [[ $GENOME == "mm10" ]]; then
  echo "Warning: mm10 is BETA (GRCm38.p4). There is no ATAQC data and blacklist (IDR peaks will not be filtered)."
  echo "Press any key to continue..."
  read -n1
else
  echo "Error: unsupported genome $GENOME"
  exit 1
fi

# URLs
if [ $GENOME == "hg19" ]; then

  CHRSZ="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.chrom.sizes"
  REF_FA="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.fa.gz"
  SEQ="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes"
  GENSZ="hs"
  BLACKLIST="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
  CHR_END_ID=22; # 22 chromosomes + chrM, chrX, chrY
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
  GENSZ="mm"
  BLACKLIST="https://personal.broadinstitute.org/anshul/projects/mouse/blacklist/mm9-blacklist.bed.gz"
  CHR_END_ID=19; # 19 chromosomes + chrM, chrX, chrY

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
  GENSZ="hs"
  BLACKLIST=""
  CHR_END_ID=22; # 22 chromosomes + chrM, chrX, chrY

elif [ $GENOME == "mm10" ]; then

  CHRSZ="http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
  REF_FA="ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M10/GRCm38.p4.genome.fa.gz"
  SEQ="http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes"
  GENSZ="mm"
  BLACKLIST=""
  CHR_END_ID=19; # 19 chromosomes + chrM, chrX, chrY
fi

# Prefix of reference genome
REF_FA_PREFIX=$(basename ${REF_FA} .gz)

# Create directories
mkdir -p ${DATA_DIR}/$GENOME
mkdir -p ${DATA_DIR}/$GENOME/seq
mkdir -p ${DATA_DIR}/$GENOME/bowtie2_index

echo "Downloading files..."
cd ${DATA_DIR}/$GENOME
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
for i in `seq 1 ${CHR_END_ID}; echo -e "M\nX\nY"`; do wget -N -c $SEQ/chr$i.fa.gz; done

echo "Extracting/processing data files..."
cd ${DATA_DIR}/$GENOME
if [ ! -f ${REF_FA_PREFIX} ]; then
  gzip -d -f -c ${REF_FA_PREFIX}.gz > ${REF_FA_PREFIX}
fi
for i in `seq 1 ${CHR_END_ID}; echo -e "M\nX\nY"`; do
  if [ ! -f seq/chr$i.fa ]; then
    gzip -f -d -c seq/chr$i.fa.gz > seq/chr$i.fa
  fi
done

echo "Building bowtie2 index..."
cd ${DATA_DIR}/$GENOME/bowtie2_index
rm -f ${REF_FA_PREFIX}
ln -s ../${REF_FA_PREFIX} ${REF_FA_PREFIX}
if [ ! -f ${REF_FA_PREFIX}.rev.1.bt2 ]; then
  source activate ${CONDA_ENV}
  bowtie2-build ${REF_FA_PREFIX} ${REF_FA_PREFIX}
fi

echo "Creating species file... (${SPECIES_FILE})"
cd ${DATA_DIR} && touch ${SPECIES_FILE}

if [[ $(grep "\[$GENOME\]" ${SPECIES_FILE} | wc -l) < 1 ]]; then
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
  echo -e "bwt2_idx\t= ${DATA_DIR}/$GENOME/bowtie2_index/${REF_FA_PREFIX}" >> ${SPECIES_FILE}
  echo -e "blacklist\t= ${BLACKLIST_PATH}" >> ${SPECIES_FILE}
  echo -e "tss_enrich\t= ${TSS_ENRICH_PATH}" >> ${SPECIES_FILE}
  echo -e "ref_fa\t= ${DATA_DIR}/$GENOME/${REF_FA_PREFIX}" >> ${SPECIES_FILE}
  echo -e "dnase\t= ${DNASE_PATH}" >> ${SPECIES_FILE}
  echo -e "prom\t= ${PROM_PATH}" >> ${SPECIES_FILE}
  echo -e "enh\t= ${ENH_PATH}" >> ${SPECIES_FILE}
  echo -e "reg2map\t= ${REG2MAP_PATH}" >> ${SPECIES_FILE}
  echo -e "roadmap_meta\t= ${ROADMAP_META_PATH}" >> ${SPECIES_FILE}
  echo "" >> ${SPECIES_FILE}
fi

# Add species file to ./default.env
sed -i "/DEF_SPECIES_FILE/c\species_file\t= ${SPECIES_FILE}\t# DEF_SPECIES_FILE: DO NOT MODIFY THIS COMMENT (install_genome_data.sh WILL NOT WORK)." "$SCRIPT_DIR/default.env"

echo "=== Genome data ($GENOME) has been successfully installed. ==="