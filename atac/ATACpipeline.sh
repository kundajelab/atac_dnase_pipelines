#!/usr/bin/env bash

set -o nounset
set -o pipefail
set -o errexit

BOWTIE_IDX=$1
READ1=$2
READ2=$3
NUMTHREADS=$4
GENOMESIZE=$5
CHROMSIZE=$6
OUTPUTDIR=$7

MAPQ_THRESH=30

mkdir -p ${OUTPUTDIR}
cd ${OUTPUTDIR}

echo "==========Trimming Adapters=========="
# trim adapters
# .fastq/.fastq.gz/.fq -> .trim.fastq
# assumes that trimAdapters outputs the R1/R2 output filenames 
# as the last 2 lines in stdout
time trimAdapters.py -a "$READ1" -b "$READ2" > >(tee trimAdapters.log) 2> >(tee trimAdapters.error.log >&2)
outputFiles=( $(tail -n 2 trimAdapters.log ) )
numOutputs=${#outputFiles[@]}
if (( numOutputs != 2 ))
then
  >&2 echo "ERROR: trimAdapters did not output the required 2 filenames"
  >&2 echo "ERROR: trimAdapters output:"
  >&2 printf "%s\n" "${outputFiles[@]}"
  >&2 echo "ERROR: --- end of output ---"
fi

trimmedR1=${outputFiles[0]}
trimmedR2=${outputFiles[1]}

echo "==========Aligning=========="
# align
time alignATAC.sh "$BOWTIE_IDX" "$trimmedR1" "$trimmedR2" "$NUMTHREADS" > >(tee alignATAC.log) 2> >(tee alignATAC.error.log >&2)
outputBAM=$(tail -n 1 alignATAC.log)

echo "==========Postprocessing=========="
# postprocess
# This output prefix assumes that the alignPostprocessPE will not overwrite
# the original BAM (i.e., will add some more suffixes other than just .bam)
postprocessBAMprefix=${outputBAM/.bam/}
time alignPostprocessPE.sh "$outputBAM" "$postprocessBAMprefix" "$MAPQ_THRESH" > >(tee alignPostprocessPE.log) 2> >(tee alignPostprocessPE.error.log >&2)
postprocessedBAM=$(tail -n 1 alignPostprocessPE.log)

echo "==========ATAC Postprocessing=========="
# ATAC specific postprocessing
readBED=${postprocessedBAM/.bam/.nonchrM.tn5.bed.gz}
time alignPostprocessATAC.sh "$postprocessedBAM" "$readBED" > >(tee alignPostprocessATAC.log) 2> >(tee alignPostprocessATAC.error.log >&2)

echo "==========Peak Calling=========="
# peak calling
time callATACpeaks.sh "$readBED" 75 "$GENOMESIZE" "$CHROMSIZE" > >(tee callATACpeaks.log) 2> >(tee callATACpeaks.error.log >&2)
