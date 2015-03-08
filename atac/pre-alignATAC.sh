#!/bin/bash

if [ -n "$1" ]; then
   SOURCE_LOG=$1
fi

outputFiles=( $(tail -n 2 ${SOURCE_LOG} ) )
numOutputs=${#outputFiles[@]}
if (( numOutputs != 2 ))
then
  >&2 echo "ERROR: trimAdapters did not output the required 2 filenames"
  >&2 echo "ERROR: trimAdapters output:"
  >&2 printf "%s\n" "${outputFiles[@]}"
  >&2 echo "ERROR: --- end of output ---"
fi

echo ${outputFiles[0]} > trimmedR1.filename
echo ${outputFiles[1]} > trimmedR2.filename
