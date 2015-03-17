#!/usr/bin/env bash

set -o pipefail
set -o errexit

if hash module 2>/dev/null; then
   if [ -n "$CLUSTER_IS_PBS" ]; then
      module add tools
      module add bowtie2/2.2.4
   else
      module add bowtie/2.2.4
   fi
   module add samtools/1.2
fi

# read from command line which files to align
if [ -n "$1" ]; then
   BOWTIE_IDX=$1 # first input is the bowtie index filename prefix (excluding trailing .X.bt2)
fi
if [ -n "$2" ]; then
   READ1=$2  # second input is read 1 fastq
fi
if [ -n "$3" ]; then
   READ2=$3  # third input is read 2 fastq
fi
if [ -n "$4" ]; then
   NUMTHREADS=$4
fi

# if these variables are the string "FROM_FILE", then load them from the
# contents of some files that should have been created previously in the pipeline
if [ "$READ1" == "FROM_FILE" ]; then
  READ1=`cat trimmedR1.filename`
fi

if [ "$READ2" == "FROM_FILE" ]; then
  READ2=`cat trimmedR2.filename`
fi

# help
if [[ -z "$READ1" ]]
then
    echo "This will align two PE fastqs"
    echo "First input is the reference"
    echo "Second input is read 1"
    echo "Third input is read 2"
    echo "Fourth input is number of threads"
    exit
fi

# if no READ2 then make it up
if [[ -z "$READ2" ]]
then
    READ2=${READ1//R1/R2}
fi

if [[ -z "$NUMTHREADS" ]]
then
    [[ -z "$NSLOTS" ]] && NSLOTS=$(nproc)
    NUMTHREADS=$NSLOTS
    >&2 echo "Using $NUMTHREADS threads for alignment"
fi

# get file type
>&2 echo "aligning" "$READ1"
type=$(basename "$READ1" | awk -F "." '{print$NF}')
inputPrefix=${READ1//_R1/}

# check if zipped
if [[ $type == 'gz' ]]
then
    # process zipped file
    outputBAM=${inputPrefix//.fastq.gz/.bam}
    outputLog=${inputPrefix//.fastq.gz/.align.log}
elif [[ $type == 'fastq' || $type == 'fq' ]]
then
    # process unzipped file
    outputBAM=${inputPrefix//.$type/.bam}
    outputLog=${inputPrefix//.$type/.align.log}
else
    >&2 echo "Unrecognized file type, accepts .fastq or .fq or .gz"
    exit 1
fi

bowtie2 -X2000 --threads "$NUMTHREADS" -x "$BOWTIE_IDX" \
    -1 <(zcat -f "$READ1") -2 <(zcat -f "$READ2") \
    2>"$outputLog" | \
    samtools view -bS - > "$outputBAM"

echo "Alignment done - output BAM file is:"
echo "$outputBAM"
