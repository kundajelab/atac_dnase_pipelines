#!/usr/bin/env bash

# this script will generate a v-plot given a sorted bam file

if hash module 2>/dev/null; then
   module add samtools/1.2
fi

if [ -n "$1" ]; then
   infile=$1 # sorted bam file
fi
if [ -n "$1" ]; then
   index=$2 # gene location bed file
fi

samtools index "$infile"

./makeVplot.py -a "$infile" -b "$index" -e 2000 -p ends -v -u
