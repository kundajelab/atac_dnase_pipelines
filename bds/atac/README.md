ATAC Seq Pipeline
===================================================


### Installation instruction

Please read this README first!
<a href="../README.md">README.md</a>


### Usage

1) Define parameters in command line argument (legacy method)
```
$ bds atac.bds [BOWTIE2_INDEX] [READ1] [READ2] [NTHREADS_BWT2] [GENOMESIZE] [CHROMSIZES] [VPLOT_INDEX] [OUTPUT_DIR]
```

2) Define parameters in command line argument.
```
$ bds atac.bds -fastq1 [READ1] -fastq2 [READ2] -bwt2_idx [BOWTIE2_INDEX] -nth_bwt2 [NTHREADS_BWT2] \
-gensz [GENOMESIZE] -chrsz [CHROMESIZES] -vplot_idx [VPLOT_INDEX]
```

3) Define parameters in configuration file.
```
$ bds atac.bds [CONF_FILE]

$ cat [CONF_FILE]
fastq1= [READ1]
fastq2= [READ2]
bwt2_idx= [BOWTIE2_INDEX]
nth_bwt2= [NTHREADS_BWT2]
gensz= [GENOMESIZE]
chrsz= [CHROMESIZES]
vplot_idx= [VPLOT_INDEX]
```

4) For Kundaje lab clusters (or using species)
Environment variables and species specific files will be automatically set in Kundaje lab clusters. Add '-kundaje_lab' at the end of the parameters.
```
$ bds atac.bds -fastq1 [READ1] -fastq2 [READ2] -kundaje_lab -species [hg19 or mm9]
```

5) For other clusters or computers
Add -mod, -addpath and -shcmd to set up enviroment variables for your jobs. This is explained later in <a href="https://github.com/kundajelab/ENCODE_chipseq_pipeline/blob/master/README_PIPELINE.md">README_PIPELINE.md</a>.


### Skipping trimming fastqs

If your fastqs are already trimmed, add the following to the commandline to skip trimming stage.
```
$ bds atac.bds ... -trimmed_fastq
```


### Processing multiple replicates IN PARALLEL

For legacy input method, # replicates is limited to 1. If you are interested in multiple replicates, don't use the legacy input method. Define fastqs with '-fastq[REP_NO]_[PAIR_NO]'.

The following example is for 3 replicates.

```
$ bds atac.bds \
-num_rep 2 \
-fastq1_1 [READ_REP1_PAIR1] \
-fastq1_2 [READ_REP1_PAIR2] \
-fastq2_1 [READ_REP2_PAIR1] \
-fastq2_2 [READ_REP2_PAIR2] \
-bwt2_idx [BOWTIE2_INDEX] -nth_bwt2 [NTHREADS_BWT2] \
-gensz [GENOMESIZE] -chrsz [CHROMESIZES] -vplot_idx [VPLOT_INDEX]
```

ATAC seq for each repliacte will go IN PARALLEL!. Consider your computation resources! # of processors taken will be :
```
max( [NTH_BWT2], [NTH_MACS2] ) x [NUM_REP]
```

If you don't want any jobs to be parallelized (each job can use multiple threads though), add '-no_par_job'.


### If you want to run pipeline on SCG3 cluster (for Kundaje lab members)

Add the following to the command line: (example for hg19):
```
$ bds -s sge atac.bds ... -kundaje_lab -species hg19_scg3 -mod "python/2.7 gnuplot/5.0"
```

Check if paths for chrsz, bwt2_idx and vplot_idx in ../chipseq/species/species_kundaje_lab.conf are valid. 
If not, 

1) modify that species file to have correct paths for them.

or

2) add the following commandline
```
-chrsz [PATH_FOR_CHRSZ_FILE] -bwt2_idx [PATH_FOR_BWT2_IDX] -vplot_idx [PATH_FOR_VPLOT_IDX]
```


### If you don't want V plot 

Add the following flag to the command line.
```
$ bds atac.bds ... -no_vplot
```


### If you don't want preseq analysis

Add the following flag to the command line.
```
$ bds atac.bds ... -no_preseq
```


### If you don't want MAC2 peak calling

Add the following flag to the command line.
```
$ bds atac.bds ... -no_peakcall
```


### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
