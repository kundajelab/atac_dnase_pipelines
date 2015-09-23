ATAC Seq Pipeline
===================================================


### Installation instruction

<a href="https://github.com/kundajelab/ENCODE_chipseq_pipeline/blob/master/README_PIPELINE.md">README_PIPELINE.md</a>

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

Add -mod, -addpath and -shcmd to set up enviroment variables for your jobs. This is explained later in this file.


### Skipping trimming fastqs

If your fastqs are already trimmed, add the following to the commandline to skip trimming stage.
```
-trimmed_fastq
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

If you don't want jobs to be parallelized (each job can use multiple threads though), add '-no_par_job'.


### If you don't want V plot 

Add the following flag to the command line.
```
-no_vplot
```

### If you don't want preseq analysis

Add the following flag to the command line.
```
-no_preseq
```


### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
