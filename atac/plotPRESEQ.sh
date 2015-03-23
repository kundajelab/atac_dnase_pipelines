#!/usr/bin/env bash

if hash module 2>/dev/null; then
   if [ -n "$CLUSTER_IS_PBS" ]; then
      module add gnuplot/5.0.0
      module add preseq/1.0.2
   fi
   module add samtools/1.2
fi

# takes unsorted bam file as input. this file should include duplicates
if [ -n "$1" ]; then
   infile=$1
fi

sortedInfile="${infile/.bam/.sort.bam}"

# preseq needs sorted input
samtools sort -Ttmp -l0 -Obam "$infile" -o "${sortedInfile}"

# run preseq
preseq lc_extrap -B -o "${sortedInfile}.preseq.dat" "${sortedInfile}" -v 2> "${sortedInfile}.preseq.log"

# plot the results
# maximum number of reads on x axis in plot (in millions)
XMAX=500

cat > preseq.gnu <<EOF
set terminal pdf 
set output '$infile.preseq.pdf'
set nokey
set style line 1 linewidth 5
set mxtics 2
set grid ytics xtics mxtics
set xrange [0:${XMAX}]
set xlabel "Number of Reads [millions]"
set ylabel "Expected Distinct Reads [millions]"
set title "PRESEQ Results for $infile"
plot '$infile.sort.bed.preseq.dat' using (\$1/1000000):(\$2/1000000) with lines linestyle 1
EOF

gnuplot preseq.gnu

