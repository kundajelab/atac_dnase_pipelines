#!/usr/bin/env bash

if hash module 2>/dev/null; then
   if [ -n "$CLUSTER_IS_PBS" ]; then
      : # texlive not installed here yet
      module add texlive/2014
   else
      module add texlive/2013
   fi
fi

# this file generates a .tex file which defines a summary report for
# this pipeline and then uses it to create report.pdf

# the first (optional) argument to the script is the data directory
if [ -n "$1" ]; then
   cd $1
fi

# the name of the folder we're working in
WD=$(pwd)
thisFolder=$(basename $WD)

# figure out what we've aligned to
MODEL=$(basename $(grep "Aligning to" alignATAC.log |awk '{print $3}'))

# insert size histogram graph file name
insertSizeHist="$(echo *.hist_graph.pdf)"

# PRESEQ graph file name
preseqGraph="$(echo *.preseq.pdf)"

# bowtie2 alignment log file name
alignLog="$(echo *.align.log)"
alignLogLastLine=$(tail -n 1 $alignLog)

# preseq log file name
preseqLog="$(echo *.preseq.log)"

# V-Plot graph file name
vPlot="$(echo *.vect.png)"

# Picard duplicate log file name
dupQCFile=$(echo *.dup.qc)
dupQCHeadings=$(sed -n 7p $dupQCFile)
dupQCHeadings=${dupQCHeadings//_/\\\\_} #moar reformatting
IFS=$'\t' read -a dupQCHeadings <<< "$dupQCHeadings"

dupQCData=$(sed -n 8p $dupQCFile)
IFS=$'\t' read -a dupQCData <<< "$dupQCData"

START=0
END=${#dupQCHeadings[@]}
DUPTABLE=""

## here we build up the tex styled table contents
for (( c=$START; c<$END; c++ ))
do
   if [ "$c" = "$((END-1))" ]; then
      DUPTABLE="${DUPTABLE}${dupQCHeadings[$c]}& ${dupQCData[$c]} \\\\ \bottomrule"
   else
      DUPTABLE="${DUPTABLE}${dupQCHeadings[$c]}& ${dupQCData[$c]} \\\\ \midrule
"
   fi
done

# library complexity data file name
bpcQCFile=$(echo *.pbc.qc)

# library complexity data file contents
declare -a bpcLabels=("TotalReadPairs" "DistinctReadPairs" "OneReadPair" "TwoReadPairs" "NRF=Distinct/Total" "PBC1=OnePair/Distinct" "PBC2=OnePair/TwoPair")
bpcQC=$(cat $bpcQCFile)
IFS=$'\t' read -a bpcData <<< "$bpcQC"

START=0
END=${#bpcLabels[@]}
BPCTABLE=""

## here we build up the tex styled table contents
for (( c=$START; c<$END; c++ ))
do
   if [ "$c" = "$((END-1))" ]; then
      BPCTABLE="${BPCTABLE}${bpcLabels[$c]}& ${bpcData[$c]} \\\\ \bottomrule"
   else
      BPCTABLE="${BPCTABLE}${bpcLabels[$c]}& ${bpcData[$c]} \\\\ \midrule
"
   fi
done

# Now we generate the tex language file that defines the report
# warning! two layers of command interpretaion in here (bash&tex),
# not for the faint of heart!
cat > "${thisFolder}.report.tex" <<EOF

\documentclass{article}
\usepackage{graphicx}
\usepackage{multicol}
\usepackage{listings}
\usepackage{grffile}
\usepackage[margin=0.5in]{geometry}
\usepackage{booktabs}
\usepackage{hyphenat}

\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\headrulewidth}{0pt}
\fancyfoot[RO, LE] {Generated on $(date)}

\begin{document}
\setlength{\columnseprule}{0.1pt}
\section{Summary for ${thisFolder//_/\\_}}
\begin{multicols}{2}
\subsection{Genome Model}
${MODEL}
\subsection{V-Plot}
From ${vPlot//_/\\_}:\\\\
\includegraphics[width=0.5\textwidth]{${vPlot}}
\subsection{preseq lc\_extrap Yield Predicion}
From ${preseqGraph//_/\\_}:\\\\
\includegraphics[width=0.5\textwidth]{${preseqGraph}}
\subsection{Library Complexity}
From ${bpcQCFile//_/\\_}:\\\\
\centerline{
\begin{tabular}{l|c}
\toprule
${BPCTABLE}
\end{tabular}
}
\subsection{Insert size histogram}
From ${insertSizeHist//_/\\_}:\\\\
\includegraphics[width=0.5\textwidth]{${insertSizeHist}}
\subsection{Picard Duplication Metrics}
From ${dupQCFile//_/\\_}:\\\\
\begin{tabular}{l|c}
\toprule
${DUPTABLE}
\end{tabular}
\subsection{bowtie2 Alignment Log}
\textbf{${alignLogLastLine//'%'/'\%'}}\\\\
From ${alignLog//_/\\_}:\\\\
\scalebox{.6}{
\lstinputlisting{${alignLog}}
}
\end{multicols}  
%\pagebreak
\end{document}
EOF

# now generate the pdf report (appears as *.report.pdf)
pdflatex "${thisFolder}.report.tex"
