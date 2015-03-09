#!/usr/bin/env bash

set -o pipefail
set -o errexit

BOWTIE_IDX=$1
READ1=$2
READ2=$3
NUMTHREADS=$4
GENOMESIZE=$5
CHROMSIZE=$6
OUTPUTDIR=$7
SCRIPTPATH=$( cd $(dirname $0) ; pwd -P ) #absolute path to directory containing this file

MAPQ_THRESH=30

mkdir -p ${OUTPUTDIR}
cd ${OUTPUTDIR}

# for cluster environment
if hash qsub 2>/dev/null; then
   echo "Cluster mode enabled"
   QSUB=true
   NEED_MEMORY=8G
   NEED_CPUS=1
   NEED_RUNTIME=20:00:00
   NEED_STACK=10M
   if hash qsig 2>/dev/null; then
      echo "PBS type cluster detected"
      PARALLEL='-l ncpus='
      MEMORY='-l mem='
      RUNTIME='-l walltime='
      STACK=''
      WAITFOR='-W afterany:'
      WD='-w'
   else
      echo "SGE type cluster detected"
      PARALLEL='-pe shm '
      MEMORY='-l h_vmem='
      RUNTIME='-l h_rt='
      STACK='-l h_stack=10M'
      WAITFOR='-hold_jid '
      WD='-wd'
   fi
fi


echo "==========Trimming Adapters=========="
# trim adapters
# .fastq/.fastq.gz/.fq -> .trim.fastq
# assumes that trimAdapters outputs the R1/R2 output filenames 
# as the last 2 lines in stdout
CMD=trimAdapters
export P1_IN="$READ1"
export P2_IN="$READ2"
if [ -n "$QSUB" ]; then
   QSUB="qsub -V ${STACK} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
${QSUB} ${SCRIPTPATH}/${CMD}.py > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)
if [ -n "$QSUB" ]; then
   JOB_ID=`head -n 1 ${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi

echo "==========Aligning=========="
# align
export SOURCE_LOG=${CMD}.log
CMD=alignATAC
if [ -n "$QSUB" ]; then
   QSUB="qsub -V -N pre-${CMD} ${WD} `pwd` ${WAITFOR}${JOB_ID} -o pre-${CMD}.log -e pre-${CMD}.error.log"
fi
${QSUB} ${SCRIPTPATH}/pre-${CMD}.sh > >(tee pre-${CMD}.log) 2> >(tee pre-${CMD}.error.log >&2)
if [ -n "$QSUB" ]; then
   JOB_ID=`head -n 1 pre-${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi

export BOWTIE_IDX=$BOWTIE_IDX
export READ1=FROM_FILE
export READ2=FROM_FILE
export NUMTHREADS=$NUMTHREADS
if [ -n "$QSUB" ]; then
   NEED_CPUS=$NUMTHREADS
   QSUB="qsub -V ${WAITFOR}${JOB_ID} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
${QSUB} ${SCRIPTPATH}/${CMD}.sh > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)
if [ -n "$QSUB" ]; then
   JOB_ID=`head -n 1 ${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi

echo "==========Postprocessing=========="
# postprocess
# This output prefix assumes that the alignPostprocessPE will not overwrite
# the original BAM (i.e., will add some more suffixes other than just .bam)
export SOURCE_LOG=${CMD}.log
CMD=alignPostprocessPE
if [ -n "$QSUB" ]; then
   QSUB="qsub -V -N pre-${CMD} ${WD} `pwd` ${WAITFOR}${JOB_ID} -o pre-${CMD}.log -e pre-${CMD}.error.log"
fi
${QSUB} ${SCRIPTPATH}/pre-${CMD}.sh > >(tee pre-${CMD}.log) 2> >(tee pre-${CMD}.error.log >&2)
if [ -n "$QSUB" ]; then
   JOB_ID=`head -n 1 pre-${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi

export RAW_BAM_FILE=FROM_FILE
export OFPREFIX=FROM_FILE
export MAPQ_THRESH=$MAPQ_THRESH
if [ -n "$QSUB" ]; then
   NEED_CPUS=1
   QSUB="qsub -V ${WAITFOR}${JOB_ID} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
${QSUB} ${SCRIPTPATH}/${CMD}.sh > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)
if [ -n "$QSUB" ]; then
   JOB_ID=`head -n 1 ${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi


echo "==========ATAC Postprocessing=========="
# ATAC specific postprocessing
export SOURCE_LOG=${CMD}.log
CMD=alignPostprocessATAC
if [ -n "$QSUB" ]; then
   QSUB="qsub -V -N pre-${CMD} ${WD} `pwd` ${WAITFOR}${JOB_ID} -o pre-${CMD}.log -e pre-${CMD}.error.log"
fi
${QSUB} ${SCRIPTPATH}/pre-${CMD}.sh > >(tee pre-${CMD}.log) 2> >(tee pre-${CMD}.error.log >&2)
if [ -n "$QSUB" ]; then
   JOB_ID=`head -n 1 pre-${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi

export input_bam=FROM_FILE
export output_file=FROM_FILE
if [ -n "$QSUB" ]; then
   NEED_CPUS=1
   QSUB="qsub -V ${WAITFOR}${JOB_ID} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
${QSUB} ${SCRIPTPATH}/${CMD}.sh > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)

echo "==========Peak Calling=========="
# peak calling
CMD=callATACpeaks

export readBed=FROM_FILE
export fragLen=75
export genomeSize=$GENOMESIZE
export chrSize=$CHROMSIZE
if [ -n "$QSUB" ]; then
   NEED_CPUS=1
   QSUB="qsub -V ${WAITFOR}${JOB_ID} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
${QSUB} ${SCRIPTPATH}/${CMD}.sh > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)
