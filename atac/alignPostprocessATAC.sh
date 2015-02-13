#!/usr/bin/env bash

# Postprocess aligned ATAC-seq BAM file
# Performs the following:
# 1. Remove reads aligning to chrM
# 2. Convert BAM to BED
# 3. Adjust ends by Tn5 insertion position
# 4. gzip output (to BED)
# 5. Generate fragment length distribution QC

set -o nounset
set -o pipefail
set -o errexit

if hash module 2>/dev/null; then
   module add bedtools/2.19.1
   module add picard-tools/1.92
   module add samtools/0.1.19
fi

input_bam=$1
output_file=$2

if [[ -z $output_file ]]
then
    output_file='/dev/stdout'
fi

# Split bam into chrM and non-chrM
nonMitoChromosomes=$(samtools view -H "${input_bam}" | \
                     grep chr | cut -f2 | sed 's/SN://g' | grep -v chrM)
nonMitoBAM=${input_bam/.bam/.nonchrM.bam}
mitoBAM=${input_bam/.bam/.chrM.bam}

samtools view -b "${input_bam}" ${nonMitoChromosomes} > "${nonMitoBAM}"
samtools view -b "${input_bam}" chrM > "${mitoBAM}"

# Process only non chrM reads
bamToBed -i "${nonMitoBAM}" | \
    adjustBedTn5.sh | \
    gzip -c > "${output_file}"

insertSizeMetrics="${PICARDROOT}/CollectInsertSizeMetrics.jar"
java -Xmx4G -jar "${insertSizeMetrics}" \
    INPUT="${nonMitoBAM}" OUTPUT="${nonMitoBAM}".hist_data.log \
    H="${nonMitoBAM}".hist_graph.pdf W=1000 STOP_AFTER=5000000
