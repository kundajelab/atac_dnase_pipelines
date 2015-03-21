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
CD="cd ${OUTPUTDIR}"
$CD

# for cluster environment
if hash qsub 2>/dev/null; then
   echo "Cluster mode enabled"
   IS_CLUSTER=true
   NEED_MEMORY=8G
   NEED_CPUS=1
   NEED_RUNTIME=5:00:00
   NEED_STACK=10M
   if hash qsig 2>/dev/null; then
      echo "PBS type cluster detected"
      export CLUSTER_IS_PBS=true
      PARALLEL='-l ncpus='
      MEMORY='-l mem='
      RUNTIME='-l walltime='
      STACK=''
      WAITFOR='-W depend=afterany:'
      WD='-w'
   else
      echo "SGE type cluster detected"
      export CLUSTER_IS_SGE=true
      PARALLEL='-pe shm '
      MEMORY='-l h_vmem='
      RUNTIME='-l h_rt='
      STACK='-l h_stack=10M'
      WAITFOR='-hold_jid '
      WD='-wd'
   fi
else
   echo "Running locally"
   COMMAND_INTERPRETER=bash
fi


echo "==========Trimming Adapters=========="
# trim adapters
# .fastq/.fastq.gz/.fq -> .trim.fastq
# assumes that trimAdapters outputs the R1/R2 output filenames 
# as the last 2 lines in stdout
CMD=trimAdapters
export P1_IN="$READ1"
export P2_IN="$READ2"
if [ -n "$IS_CLUSTER" ]; then
   COMMAND_INTERPRETER="qsub -V ${STACK} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
echo "$CD; ${SCRIPTPATH}/${CMD}.py" | $COMMAND_INTERPRETER > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)
if [ -n "$IS_CLUSTER" ]; then
   JOB_ID=`head -n 1 ${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi

echo "==========Aligning=========="
# align
export SOURCE_LOG=${CMD}.log
CMD=alignATAC
export BOWTIE_IDX=$BOWTIE_IDX
export READ1=FROM_FILE
export READ2=FROM_FILE
export NUMTHREADS=$NUMTHREADS
if [ -n "$IS_CLUSTER" ]; then
   NEED_CPUS=$NUMTHREADS
   NEED_RUNTIME=10:00:00
   if [ -n "$CLUSTER_IS_SGE" ]; then
      NEED_MEMORY=2G
   fi
   COMMAND_INTERPRETER="qsub -V ${WAITFOR}${JOB_ID} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
echo "$CD; ${SCRIPTPATH}/pre-${CMD}.sh; ${SCRIPTPATH}/${CMD}.sh" | $COMMAND_INTERPRETER > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)
if [ -n "$IS_CLUSTER" ]; then
   JOB_ID=`head -n 1 ${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi

echo "==========Postprocessing=========="
# postprocess
# This output prefix assumes that the alignPostprocessPE will not overwrite
# the original BAM (i.e., will add some more suffixes other than just .bam)
export SOURCE_LOG=${CMD}.log
CMD=alignPostprocessPE
export RAW_BAM_FILE=FROM_FILE
export OFPREFIX=FROM_FILE
export MAPQ_THRESH=$MAPQ_THRESH
if [ -n "$IS_CLUSTER" ]; then
   NEED_CPUS=1
   if [ -n "$CLUSTER_IS_SGE" ]; then
      NEED_MEMORY=8G
   fi
   COMMAND_INTERPRETER="qsub -V ${WAITFOR}${JOB_ID} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
echo "$CD; ${SCRIPTPATH}/pre-${CMD}.sh; ${SCRIPTPATH}/${CMD}.sh" | $COMMAND_INTERPRETER > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)
if [ -n "$IS_CLUSTER" ]; then
   JOB_ID=`head -n 1 ${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi


echo "==========ATAC Postprocessing=========="
# ATAC specific postprocessing
export SOURCE_LOG=${CMD}.log
CMD=alignPostprocessATAC
export input_bam=FROM_FILE
export output_file=FROM_FILE
if [ -n "$IS_CLUSTER" ]; then
   COMMAND_INTERPRETER="qsub -V ${WAITFOR}${JOB_ID} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
echo "$CD; ${SCRIPTPATH}/pre-${CMD}.sh; ${SCRIPTPATH}/${CMD}.sh" | $COMMAND_INTERPRETER > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)
if [ -n "$IS_CLUSTER" ]; then
   JOB_ID=`head -n 1 ${CMD}.log | awk 'match($0,/[0-9]+/){print substr($0, RSTART, RLENGTH)}'`
fi


echo "==========Peak Calling=========="
# peak calling
CMD=callATACpeaks
export readBed=FROM_FILE
export fragLen=75
export genomeSize=$GENOMESIZE
export chrSize=$CHROMSIZE
if [ -n "$IS_CLUSTER" ]; then
   COMMAND_INTERPRETER="qsub -V ${WAITFOR}${JOB_ID} -N ${CMD} ${WD} `pwd` ${MEMORY}${NEED_MEMORY} ${PARALLEL}${NEED_CPUS} ${RUNTIME}${NEED_RUNTIME} -o ${CMD}.log -e ${CMD}.error.log"
fi
echo "$CD; ${SCRIPTPATH}/${CMD}.sh" | $COMMAND_INTERPRETER > >(tee ${CMD}.log) 2> >(tee ${CMD}.error.log >&2)
