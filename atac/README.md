ATAC-seq pipeline
=================

End-to-end ATAC-seq pipeline from FASTQ to peaks.

Input: FASTQs (.fastq or .fastq.gz)

Performs the following steps:
* Adapter trimming (using Jason's script)
* Alignment (using Bowtie2)
* Postprocessing / QC (generates insert size distribution)
* Peak calling using MACS2 (using the shifting trick)

Usage
=====

./alignATAC BOWTIE_IDX READ1 READ2 NUMTHREADS GENOMESIZE CHROMSIZE OUTPUTDIR

BOWTIE_IDX: path to BOWTIE2 index
READ1: read 1 fastq(.gz)
READ2: read 2 fastq(.gz)
NUMTHREADS: number of threads for bowtie2
GENOMESIZE: hs / mm (this is for MACS2)
CHROMSIZE: chromosome size file (may be found in $GENOMESIZEDIR)
OUTPUTDIR: directory to put MACS2 peak calls / signal tracks
