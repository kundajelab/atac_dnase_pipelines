#!/bin/bash

if [ -n "$1" ]; then
   SOURCE_LOG=$1
fi

outputBAM=$(tail -n 1 ${SOURCE_LOG})
postprocessBAMprefix=${outputBAM/.bam/}

echo ${outputBAM} > outputBAM.filename
echo ${postprocessBAMprefix} > postprocessBAMprefix.filename
