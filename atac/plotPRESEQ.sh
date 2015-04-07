#!/usr/bin/env bash

if hash module 2>/dev/null; then
   if [ -n "$CLUSTER_IS_PBS" ]; then
      module add gnuplot/5.0.0
   else
      module add gnuplot/5.0
   fi
   module add preseq/1.0.2
fi

# takes sorted bam file as input. this file should include duplicates
if [ -n "$1" ]; then
   sortedInfile=$1
fi

# run preseq
preseqData="${sortedInfile}.preseq.dat"
preseq lc_extrap -P -B -o "${preseqData}" "${sortedInfile}" -v 2> "${sortedInfile}.preseq.log"

# plot the results
# maximum number of reads on x axis in plot (in millions)
XMAX=500

titleName=$(basename $sortedInfile) 
titleName=${titleName//_/\\_ } #escape underscores
cat > preseq.gnu <<EOF
set terminal pdfcairo
set output '$sortedInfile.preseq.pdf'
set key box bottom right
set style line 1 linewidth 5
set style line 2 linewidth 1
set mxtics 2
set grid ytics xtics mxtics
set xrange [0:${XMAX}]
set xlabel "Number of Reads [millions]"
set ylabel "Expected Distinct Reads [millions]"
set title "PRESEQ Results for $titleName"
plot '${preseqData}' using (\$1/1000000):(\$2/1000000) with lines linestyle 1 notitle, \
                  '' using (\$1/1000000):(\$3/1000000) with lines linestyle 2 notitle, \
                  '' using (\$1/1000000):(\$4/1000000) with lines linestyle 2 title '+/- 95% confidence interval'
EOF

gnuplot preseq.gnu
