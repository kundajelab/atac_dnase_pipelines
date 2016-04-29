ATAC Seq Pipeline
===================================================


### Installation instruction (for SCG3 and all others)

Install java (jdk >= 1.7 or jre >= 1.7) and the latest git on your system. 

Install Anaconda Python3 (or Miniconda3) on your system. If you already have Anaconda Python 3, skip it. Open a new terminal after installation.

Install BigDataScript v0.9999 on your system.

See details <a href="https://github.com/kundajelab/TF_chipseq_pipeline/blob/master/README_PIPELINE.md" target=_blank>here</a>

ATAC Seq pipeline is dependent on two external git repos and each has the following directory:
```
ataqc/ 		(public, https://bitbucket.org/csfoo/ataqc)
chipseq/	(public, https://github.com/kundajelab/TF_chipseq_pipeline)
```

To FULLY clone the repo:
```
$ git clone https://github.com/kundajelab/bds_atac --recursive
```

To FULLY pull the repo:
```
$ git pull --recurse-submodules
```

Install software dependencies automatically (DO NOT run this on kundaje clusters). It will create two conda environments (aquas_chipseq and aquas_chipseq_py3) in Miniconda3.
```
$ ./install_dependencies.sh
```

Replace BDS's default `bds.config` with a correct one:
```
$ cp bds.config $HOME/.bds
```


### Installation instruction (for Kundaje lab clusters)

BDS and all dependencies have already been installed on Kundaje lab servers. Do not run `install_dependencies.sh` on these servers. Get the latest version of chipseq pipelines. Don't forget to move bds.config to BigDataScript (BDS) directory
```
$ git clone https://github.com/kundajelab/bds_atac --recursive
$ cd bds_atac
$ mkdir -p $HOME/.bds
$ cp chipseq/bds.config $HOME/.bds/
```


### Genome data files for SCG3 and Kundaje lab servers

For Kundaje lab servers (mitra, nandi, durga, kali, vayu, amold and wotan) and SCG3 (carmack*, crick*, scg3*), the pipeline automatically determines the type of servers and set shell environments and species database.
```
$ bds bds_atac.bds ... -species [SPECIES; hg19, mm9, ... ]
```


### Usage


1) Define parameters in command line argument.
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

For ATAQC, you need to define the following parameters. parameters `-preseq` and `-vplot` will be ignored since they are already included in ATAQC. See help (`$ bds atac.bds`) for description of all parameters. Even though you don't use a species file `-species_file`, you need to specify a species name for ATAQC.
```
-species [hg19, mm9 or ...] -tss_enrich [] -ref_fa [] -blacklist [] -dnase [] -prom [] -enh [] -reg2map [] -roadmap_meta []
```

If you want to just align your data (no peak calling or further steps like IDR).
```
-align
```

If you don't want ATAQC, add the following to command line. 
```
-no_ataqc 
```

If you have just one replicate (PE), define fastqs with `-fastq[PAIR_NO]`.
```
-fastq1 [READ_PAIR1] -fastq2 [READ_PAIR2] \
```

For multiple replicates (PE), define fastqs with `-fastq[REP_NO]_[PAIR_NO]`. Add -fastq[]_[] for each replicate and pair to the command line:replicates.
```
-fastq1_1 [READ_REP1_PAIR1] -fastq1_2 [READ_REP1_PAIR2] -fastq2_1 [READ_REP2_PAIR1] -fastq2_2 [READ_REP2_PAIR2] ...
```

For multiple replicates (SE), define fastqs with `-fastq[REP_NO]`:
```
-se -fastq1 [READ_REP1] -fastq2 [READ_REP2] ...
```

You can also start from bam files. There are two kinds of bam files (raw or deduped) and you need to explicitly choose between raw bam (bam) and deduped one (nodup_bam) with `-input [BAM_TYPE]`. Don't forget to add `-se` if they are not paired end (PE).

For raw bams,
```
-bam1 [RAW_BAM_REP1] -bam2 [RWA_BAM_REP1] ...
```

For deduped (filtered) bams, preseq analysis and TSS enrichment plot will not be available since they need sorted raw bam.
```
-filt_bam1 [NODUP_BAM_REP1] -filt_bam2 [NODUP_BAM_REP1] ...
```

For tagaligns (non-tn5-shifted), preseq analysis and TSS enrichment plot will not be available since they need sorted raw bam.
```
-tag1 [TAGALIGN_REP1] -tag2 [TAGALIGN_REP2] ...
```

To subsample beds (tagaligns) add the following to the command line. This is different from subsampling for cross-corr. analysis. Peaks will be called with subsampled tagaligns.
```
-subsample [NO_READS_TO_SUBSAMPLE]
```

To change # of lines to subsample for cross-corr. analysis.
```
-nreads [NO_READS_TO_SUBSAMPLE]
```

To disable pseudo replicate generation. By default, IDR will be done for true replicates and pseudo replicates, but if you have `-true_rep` in the command line, you will also get IDR on true replicates only. IDR on a single replicate and naive overlapped peak is not avaiable when this flag is on:
```
-true_rep
```

IDR analysis is included in the pipeline by default. For better IDR QC, add path to blacklist idr (for hg19, http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz).
```
-blacklist_idr [BLACKLIST_IDR]
```

If you don't want IDR analysis on peaks (two replicates are needed) add the following:
```
-no_idr
```

To change resource settings (# of processor, max memory and walltime) for bowtie2, add the following to command line, please note that memory is PER CPU:
```
-nth_bwt2 [NTHREADS_BWT2] -mem_bwt2 [MAX_MEMORY_PER_THREAD_BWT2; e.g. 20G] -wt_bwt2 [WALLTIME_BWT2; e.g. 20h]
```

For MACS2 peak calling:
```
-nth_macs2 [NTHREADS_MACS2] -mem_macs2 [MEMORY_MACS2; e.g. 20G] -wt_macs2 [WALLTIME_MACS2; e.g. 20h]
```

For Kundaje lab cluster and SCG3, skip all genome specific parameters (like bwt2_idx, chrsz, ... ) and just specify species.
```
$ bds atac.bds -fastq1 [READ1] -fastq2 [READ2] -species [hg19 or mm9]
```

For other clusters, add -mod, -addpath, -shcmd, -conda_env to set up enviroment variables for your jobs or make an environment file for your system. See details <a href="https://github.com/kundajelab/ENCODE_chipseq_pipeline/blob/master/README_PIPELINE.md">here</a>.

To list all parameters and default values for them,
```
$ bds atac.bds
```


2) Define parameters in configuration file.
Key names in a configruation file are identical to parameter names on command line. 
```
$ bds atac.bds [CONF_FILE]

or 

$ bds atac.bds -c [CONF_FILE]

$ cat [CONF_FILE]
fastq1= [READ1]
fastq2= [READ2]
...
```


### Species file and Environment file

See details <a href="https://github.com/kundajelab/TF_chipseq_pipeline/blob/master/README_PIPELINE.md" target=_blank>here</a>



### Parallelization level

ATAC seq for each repliacte will go IN PARALLEL!. Consider your computation resources! # of processors taken will be :
```
max( [NTH_BWT2], [NTH_MACS2], [NTH_SPP] ) x [NUM_REP]
```

For completely serialized jobs:
```
-no_par
```

You can also set up the level of parallelization for the pipeline.
```
-par_lvl [PAR_LEVEL; 0-7]
```
0: no parallel jobs (equivalent to `-no_par`, all subtasks for each replicate will also be serialized)
1: no replicates/controls in parallel (subtasks for each replicate can be parallelized)
2: 2 replicates/controls in parallel
3: 2 replicates/controls and 2 peak-callings in parallel (default)
4: 4 replicates/controls and 2 peak-callings in parallel
5: 4 replicates/controls and 4 peak-callings in parallel
6: customized
7: unlimited

For customized parallelization:
```
-par_lvl 6 -reps_in_par [NO_REP_IN_PAR] -peaks_in_par [NO_PEAKCALLING_IN_PAR]
```

See details <a href="https://github.com/kundajelab/TF_chipseq_pipeline/blob/master/README_PIPELINE.md" target=_blank>here</a>



### Requirements 

For python2 (python 2.x >= 2.7), <a href="https://github.com/kundajelab/bds_atac/blob/master/requirements.txt" target=_blank>here</a>

For python3, <a href="https://github.com/kundajelab/bds_atac/blob/master/requirements_py3.txt" target=_blank>here</a>

For R-2.x, <a href="https://github.com/kundajelab/bds_atac/blob/master/requirements_r2.txt" target=_blank>here</a>



### Troubleshooting

1) [main_samview] random alignment retrieval only works for indexed BAM or CRAM files.

If your pipeline starts from BAM files, make sure that bam index (.bam.bai) exists together with BAM file. If not, build index with `samtools index [YOUR_BAM_FILE]`. BAM and BAI should be in the same directory.


2) pysam backward compatibility issue

ATAQC currently does not work with pysam >= 0.9. Lower it to 0.8.3.
```
Traceback (most recent call last):
  File "/users/leepc12/code/bds_atac/ataqc/run_ataqc.py", line 1303, in <module>
    main()
  File "/users/leepc12/code/bds_atac/ataqc/run_ataqc.py", line 1120, in main
    chr_m_reads, fraction_chr_m = get_chr_m(COORDSORT_BAM)
  File "/users/leepc12/code/bds_atac/ataqc/run_ataqc.py", line 160, in get_chr_m
    tot_reads += int(chrom_stats[2])
IndexError: list index out of range
```

3) samtools ncurses bug

Prepend a directory for `libncurses.so.5` to `LD_LIBRARY_PATH`. See `install_dependencies.sh` for solution.
```
samtools: symbol lookup error: /lib/x86_64-linux-gnu/libncurses.so.5: undefined symbol: _nc_putchar
```

4) Error: could not find environment: bds_atac

Unload any Anaconda Python modules. Remove locally installed Anaconda Python from your `$PATH`



### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Chuan Sheng Foo - PhD Student, Computer Science Dept., Stanford University
* Daniel Kim - MD/PhD Student, Biomedical Informatics Program, Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
