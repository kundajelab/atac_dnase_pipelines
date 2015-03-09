#!/bin/bash

if [ -n "$1" ]; then
   SOURCE_LOG=$1
fi

postprocessedBAM=$(tail -n 1 ${SOURCE_LOG})
readBED=${postprocessedBAM/.bam/.nonchrM.tn5.bed.gz}

echo ${postprocessedBAM} > postprocessedBAM.filename
echo ${readBED} > readBED.filename
