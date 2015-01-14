ATAC-seq pipeline
=================

End-to-end ATAC-seq pipeline from FASTQ to peaks.

Input: FASTQs (.fastq or .fastq.gz)

Performs the following steps:
* Adapter trimming (using Jason's script)
* Alignment (using Bowtie2)
* Postprocessing / QC (generates insert size distribution)
* Peak calling using MACS2 (using the shifting trick)
