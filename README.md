ATAC Seq Pipeline
===================================================



### Cloning/pulling the pipeline repo

ATAC Seq pipeline is dependent on two external git repos and each has the following directory:
```
ataqc/ 		(public, https://bitbucket.org/csfoo/ataqc)
chipseq/	(public, https://github.com/kundajelab/TF_chipseq_pipeline)
```

To fully clone the repo:
```
$ git clone https://github.com/kundajelab/bds_atac --recursive
```

To fully pull the repo:
```
$ git pull --recurse-submodules
```


### Installing BigDataScript (BDS)

Get BigDataScript v0.9999:
```
$ git clone https://github.com/pcingola/BigDataScript
$ cd BigDataScript
$ git checkout tags/v0.9999
$ cp distro/bds_Linux.tgz $HOME
$ cd $HOME
$ tar zxvf bds_Linux.tgz
```

Add `$HOME/.bds/` to your `$PATH`, then replace BDS's default bds.config with a correct one:
```
$ cp /path/to/bds_atac/chipseq/bds.config $HOME/.bds/
```

If Java memory error occurs, add the following to your `$HOME/.bashrc`:
```
export _JAVA_OPTIONS="-Xms256M -Xmx512M -XX:ParallelGCThreads=1"
export MAX_JAVA_MEM="8G"
export MALLOC_ARENA_MAX=4
```





### Usage

1) Define parameters in command line argument (legacy method)
This input method does not support multiple replicates and always generate V plot and perform preseq analysis.
```
$ bds atac.bds [BOWTIE2_INDEX] [READ1] [READ2] [NTHREADS_BWT2] [GENOMESIZE; hs for human, mm for mouse] [CHROMSIZES_FILE] [VPLOT_INDEX] [OUTPUT_DIR]
```

2) Define parameters in command line argument.
For general use, use the following command line:
```
$ bds atac.bds -fastq1 [READ1] -fastq2 [READ2] -bwt2_idx [BOWTIE2_INDEX] \
-gensz [GENOMESIZE; hs for human, mm for mouse] -chrsz [CHR_SIZES_FILE]
```

If your fastqs are already trimmed, add the following to the command line to skip trimming stage.
```
-trimmed_fastq
```

If your data are single ended, add the following to the command line.
```
-se
```

For V plot generation, add the following to command line:
```
-vplot -vplot_idx [VPLOT_INDEX] 
```

For preseq analysis, add the following to command line:
```
-preseq
```

For advanced ATAQC (only for PE dataset), add the following to command line, parameters `-preseq` and `-vplot` will be ignored since they are already included in ATAQC. You will need to get read permission to the ataqc repo (https://bitbucket.org/csfoo/ataqc).
```
-ataqc -vplot_idx [VPLOT_INDEX]
```

If you have just one replicate (PE), define fastqs with `-fastq[PAIR_NO]`.
```
-fastq1 [READ_PAIR1] -fastq2 [READ_PAIR2] \
```

For multiple replicates (PE), define fastqs with `-fastq[REP_NO]_[PAIR_NO]`. Add -fastq[]_[] for each replicate and pair to the command line:replicates.
```
-fastq1_1 [READ_REP1_PAIR1] -fastq1_2 [READ_REP1_PAIR2] \
-fastq2_1 [READ_REP2_PAIR1] -fastq2_2 [READ_REP2_PAIR2] \
...
```

For multiple replicates (SE), define fastqs with `-fastq[REP_NO]`:
```
-se \
-fastq1 [READ_REP1] \
-fastq2 [READ_REP2] \
...
```


You can also start from bam files. There are two kinds of bam files (raw or deduped) and you need to explicitly choose between raw bam (bam) and deduped one (nodup_bam) with `-input [BAM_TYPE]`.

For raw bams,
```
-bam1 [RAW_BAM_REP1] -bam2 [RWA_BAM_REP1] ...
```

For deduped (filtered) bams, preseq analysis and v plot will not be available since they need sorted raw bam.
```
-filt_bam1 [NODUP_BAM_REP1] -filt_bam2 [NODUP_BAM_REP1] ...
```

To subsample beds (tagaligns) add the following to the command line, you can skip the second parameter (-nreads, default is 15000000):
```
-subsample -nreads [NO_READS_TO_SUBSAMPLE]
```

To generate pseduro replicates and call peaks on them:
```
-pseudorep
```

For IDR analysis on peaks (two replicates are needed):
```
-idr
```

To change resource settings (# of processor, max memory and walltime) for bowtie2, add the following to command line:
```
-nth_bwt2 [NTHREADS_BWT2] -mem_bwt2 [MEMORY_BWT2; e.g. 20G] -wt_bwt2 [WALLTIME_BWT2; e.g. 20h]
```

For MACS2 peak calling:
```
-nth_macs2 [NTHREADS_MACS2] -mem_macs2 [MEMORY_MACS2; e.g. 20G] -wt_macs2 [WALLTIME_MACS@; e.g. 20h]
```


By default, IDR will be done for true replicates, but if you have `-pseudorep` in the command line, you will also get IDR on pseudo replicates and pooled pseudo replicates.


For Kundaje lab cluster and SCG3, skip parameters (bwt2_idx, chrsz, gensz and vplot_idx) and just specify species.
```
$ bds atac.bds -fastq1 [READ1] -fastq2 [READ2] -species [hg19 or mm9]
```

For other clusters, add -mod, -addpath and -shcmd to set up enviroment variables for your jobs. This is explained in <a href="https://github.com/kundajelab/ENCODE_chipseq_pipeline/blob/master/README_PIPELINE.md">README_PIPELINE.md</a>.

To list all parameters and default values for them,
```
$ bds atac.bds
```



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
-no_par
```

### Requirements (python 2.x >= 2.7)

Python modules:
- numpy
- matplotlib
- pysam
- python-Levenshtein
- pybedtools
- trimgalore



### Troubleshooting

1) [main_samview] random alignment retrieval only works for indexed BAM or CRAM files.

If your pipeline starts from BAM files, make sure that bam index (.bam.bai) exists together with BAM file. If not, build index with `samtools index [YOUR_BAM_FILE]`. BAM and BAI should be in the same directory.



### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Chuan Sheng Foo - PhD Student, Computer Science Dept., Stanford University
* Daniel Kim - MD/PhD Student, Biomedical Informatics Program, Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
