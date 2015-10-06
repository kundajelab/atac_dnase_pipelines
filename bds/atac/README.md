ATAC Seq Pipeline
===================================================


### Installation instruction

Please read this README first!
<a href="../README.md">README.md</a>


### Usage

1) Define parameters in command line argument (legacy method)
This input method does not support multiple replicates and always generate V plot and perform preseq analysis.
```
$ bds atac.bds [BOWTIE2_INDEX] [READ1] [READ2] [NTHREADS_BWT2] [GENOMESIZE] [CHROMSIZES] [VPLOT_INDEX] [OUTPUT_DIR]
```

2) Define parameters in command line argument.
For general use, use the following command line:
```
$ bds atac.bds -fastq1 [READ1] -fastq2 [READ2] -bwt2_idx [BOWTIE2_INDEX] -nth_bwt2 [NTHREADS_BWT2] \
-gensz [GENOMESIZE] -chrsz [CHROMESIZES]
```

If your fastqs are already trimmed, add the following to the command line to skip trimming stage.
```
-trimmed_fastq
```

For V plot generation, add the following to command line:
```
-vplot -vplot_idx [VPLOT_INDEX] 
```

For preseq analysis, add the following to command line:
```
-preseq
```

For multiple replicates, define fastqs with '-fastq[REP_NO]_[PAIR_NO]'. Add -num_rep and -fastq[]_[] for each replicate and pair to the command line. The following example is for 2 replicates.
```
-num_rep 2 \
-fastq1_1 [READ_REP1_PAIR1] -fastq1_2 [READ_REP1_PAIR2] \
-fastq2_1 [READ_REP2_PAIR1] -fastq2_2 [READ_REP2_PAIR2]
```

To change resource settings (# of processor, max memory and walltime) for bowtie2, add the following to command line:
```
-nth_bwt2 [NTHREADS_BWT2] -mem_bwt2 [MEMORY_BWT2; e.g. 20G] -wt_bwt2 [WALLTIME_BWT2; e.g. 20h]
```

For Kundaje lab cluster and SCG3, skip parameters (bwt2_idx, chrsz, gensz and vplot_idx) and just specify species.
```
$ bds atac.bds -fastq1 [READ1] -fastq2 [READ2] -species [hg19 or mm9]
```

For other clusters, add -mod, -addpath and -shcmd to set up enviroment variables for your jobs. This is explained in <a href="https://github.com/kundajelab/ENCODE_chipseq_pipeline/blob/master/README_PIPELINE.md">README_PIPELINE.md</a>.



3) Define parameters in configuration file.
Key names in a configruation file are identical to parameter names on command line. 
```
$ bds atac.bds [CONF_FILE]

$ cat [CONF_FILE]
fastq1= [READ1]
fastq2= [READ2]
...
```


### Parallelism in atac pipeline

ATAC seq for each repliacte will go IN PARALLEL!. Consider your computation resources! # of processors taken will be :
```
max( [NTH_BWT2], [NTH_MACS2] ) x [NUM_REP]
```

(Not recommended) If you don't want any jobs to be parallelized (each job can still use multiple threads though), add the following to command line (this option is for computers with limited resource):
```
-no_par_job
```


### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
