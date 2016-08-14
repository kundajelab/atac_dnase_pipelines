ATAC-Seq / DNase-Seq Pipeline
===================================================


# Installation

<a href="https://github.com/kundajelab/bds_atac/blob/master/INSTALL.md">General</a>

<a href="https://github.com/kundajelab/bds_atac/blob/master/INSTALL_Kundaje.md">Kundaje lab</a>

<a href="https://github.com/kundajelab/bds_atac/blob/master/INSTALL_SCG_Sherlock.md">SCG3/4 and Stanford Sherlock clusters</a>


# Usage

We recommend using BASH to run pipelines.

For general use, use the following command line: (for PE data set)
```
$ bds atac.bds -species [SPECIES; hg19, mm9, ... ] -nth [NUM_THREADS] -fastq1_1 [READ1] -fastq1_2 [READ2]
```

For DNase-seq: (it's <b>NOT `-dnase`</b>!)
```
-dnase_seq
```

The pipeline can automatically detect adapters if you remove `-adapter` from your command line. (<b>You need to have an access to the git repo https://github.com/nboley/GGR_code</b>)

To use old adapter trimmers (`trim_galore` for SE and `trimAdapter.py` for PE):
```
-old_trimmer
```

If your fastqs are already trimmed, add the following to the command line to skip trimming stage.
```
-trimmed_fastq
```
<b>IMPORTANT!</b> If your data set is <b>SINGLE ENDED</b> add the following to the command line, otherwise the pipeline works for PE by default.
```
-se 
```
You can also individually specify endedness for each replicate.
```
-se[REPLICATE_ID] 	# for exp. replicates, 
```
If you want to just align your data (no peak calling or further steps like IDR).
```
-align
```
If you don't want ATAQC, add the following to command line. 
```
-no_ataqc 
```
If you have just one replicate (PE), define fastqs with `-fastq[REP_ID]_[PAIR_ID]`.
```
-fastq1_1 [READ_PAIR1] -fastq1_2 [READ_PAIR2] \
```
For multiple replicates (PE), define fastqs with `-fastq[REP_ID]_[PAIR_ID]`. Add `-fastq[]_[]` for each replicate and pair to the command line:replicates.
```
-fastq1_1 [READ_REP1_PAIR1] -fastq1_2 [READ_REP1_PAIR2] -fastq2_1 [READ_REP2_PAIR1] -fastq2_2 [READ_REP2_PAIR2] ...
```
For multiple replicates (SE), define fastqs with `-fastq[REP_ID]`:
```
-se -fastq1 [READ_REP1] -fastq2 [READ_REP2] ...
```
You can also specify an adapter to be trimmed for each fastq. Example: `-adapter1_1 [ADAPTER1_1] -adapter1_2 [ADAPTER1_2] ...` for PE or `-adapter1 [ADAPTER1] -adapter2 [ADAPTER2] ...`. Define adapters just like how you defined your fastqs.

You can start from bam files. There are two kinds of bam files (raw or deduped) and you need to explicitly choose between raw bam (bam) and deduped one (filt_bam) with `-input [BAM_TYPE]`. Don't forget to add `-se` if they are not paired end (PE). For raw bams,
```
-bam1 [RAW_BAM_REP1] -bam2 [RWA_BAM_REP1] ...
```
For deduped (filtered) bams:
```
-filt_bam1 [FILT_BAM_REP1] -filt_bam2 [FILT_BAM_REP1] ...
```
For tagaligns (non-tn5-shifted):
```
-tag1 [TAGALIGN_REP1] -tag2 [TAGALIGN_REP2] ...
```
You can also mix up any data types.
```
-bam1 [RAW_BAM_REP1] -tag2 [TAGALIGN_REP2] ...
```
You can also define endedness (SE/PE) of each replicate. (example: SE for replicate 1, PE for replicate 2, SE for replicate 3):
```
-se1 -pe2 -se3
```
To subsample beds (tagaligns) add the following to the command line. This is different from subsampling for cross-corr. analysis. Peaks will be called with subsampled tagaligns.
```
-subsample [NO_READS_TO_SUBSAMPLE]
```
To change # of lines to subsample for cross-corr. analysis. This will not affect tasks downstream (peak calling and IDR).
```
-nreads [NO_READS_TO_SUBSAMPLE]
```
To disable pseudo replicate generation, add the following. By default, peak calling and IDR will be done for true replicates and pseudo replicates, but if you have `-true_rep` in the command line, you will also get IDR on true replicates only.
```
-true_rep
```
IDR analysis is included in the pipeline by default. If there are more than two replicates, IDR will be done for every pair of replicates. If you don't want IDR, add the following:
```
-no_idr
```
To list all parameters and default values for them,
```
$ bds atac.bds
```

You can also define parameters in a configuration file. Key names in a configruation file are identical to parameter names on command line. 
```
$ bds atac.bds [CONF_FILE]

$ cat [CONF_FILE]
species = [SPECIES; hg19, mm9, ...]
nth 	= [NUM_THREADS]
fastq1_1= [READ1]
fastq1_2= [READ2]
...
```

## Parallelization

For completely serialized jobs, add `-no_par` to the command line. Individual tasks can still go multi-threaded.

**IMPORTANT!** You can set up a limit for total number of threads with `-nth [MAX_TOTAL_NO_THREADS]`. Total number of threads used by a pipeline will not exceed this limit.

Default `-nth` for each cluster is defined on `./default.env` (e.g. 16 on SCG and 8 on Kundaje lab cluster)

The pipeline automatically distributes `[MAX_TOTAL_NO_THREADS]` threads for jobs according to corresponding input file sizes. For example of two fastqs (1GB and 2GB) with `-nth 6`, 2 and 4 threads are allocated for aligning 1GB and 2GB fastqs, respectively. The same policy applies to other multi-threaded tasks like deduping and peak calling.

However, all multi-threaded tasks (like bwa, bowtie2, spp and macs2) still have their own max. memory (`-mem_APPNAME [MEM_APP]`) and walltime (`-wt_APPNAME [WALLTIME_APP]`) settings. Max. memory is **NOT PER CPU**. For example on Kundaje lab cluster (with SGE flag activated `bds -s sge bds_atac.bds ...`) or on SCG4, the actual shell command submitted by BDS for each task is like the following:

```
qsub -V -pe shm [NTH_ALLOCATED_FOR_APP] -h_vmem=[MEM_APP]/[NTH_ALLOCATED_FOR_APP] -h_rt=[WALLTIME_APP] -s_rt=[WALLTIME_APP] ...
```

This ensures that total memory reserved for a cluster job equals to `[MEM_APP]`. The same policy applies to SLURM.

## Managing multiple pipelines

`./utils/bds_scr` is a BASH script to create a detached screen for a BDS script and redirect stdout/stderr to a log file `[LOG_FILE_NAME]`. If a log file already exists, stdout/stderr will be appended to it.

Monitor the pipeline with `tail -f [LOG_FILE_NAME]`. The only difference between `bds_scr` and `bds` is that you have `[SCR_NAME] [LOG_FILE_NAME]` between `bds` command and its parameters (or a BDS script name).

You can skip `[LOG_FILE_NAME]` then a log file `[SCR_NAME].log` will be generated on the working directory.

You can also add any BDS parameters (like `-dryRun`, `-d` and `-s`). The following example is for running a pipeline on Sun Grid Engine.

```
$ bds_scr [SCR_NAME] [LOG_FILE_NAME] bds_atac.bds ...
$ bds_scr [SCR_NAME] bds_atac.bds ...
$ bds_scr [SCR_NAME] -s sge bds_atac.bds ...
```

Once the pipeline run is done, the screen will be automatically closed. To kill a pipeline manually while it's running, use `./utils/kill_scr` or `screen -X quit`:

```
$ screen -X -S [SCR_NAME] quit
$ kill_scr [SCR_NAME]
```

## Useful HTML reports

There are two kinds of HTML reports provided by the pipeline.

* BigDataScript HTML report for debugging: Located at the working folder with name bds_atac_[TIMESTAMP]_report.html. This report is automatically generated by BigDataScript and is useful for debugging since it shows summary, timeline, Stdout and Stderr for each job.

* ATAC-Seq pipeline report for QC and result: The pipeline automatically generate a nice HTML report (Report.html) on its output directory (specified with -out_dir or just './out'). It summarizes files and directory structure, includes QC reports and show a workflow diagram and genome browser tracks for peaks and signals (bigwigs for pValue and fold change). Move your output directory to a web directory (for example, /var/www/somewhere) or make a softlink of it to a web directory. For genome browser tracks, specify your web directory root for your output  While keeping its structure. Make sure that you have bgzip and tabix installed on your system. Add the following to the command line:

      -url_base http://your/url/to/output -title [PREFIX_FOR_YOUR_REPORT]

## Coding with BDS

* [Using modules in AQUAS pipeline](https://github.com/kundajelab/TF_chipseq_pipeline/blob/master/README_CODE.md)

* [BigDataScript github repo](https://github.com/pcingola/BigDataScript)

* [BigDataScript documentation](http://pcingola.github.io/BigDataScript/bigDataScript_manual.html)

# Requirements 

* For python2 (python 2.x >= 2.7) and R-3.x, [requirements.txt](requirements.txt)

* For python3, [requirements_py3.txt](requirements_py3.txt)

# Troubleshooting

See [more troubleshootings](https://github.com/kundajelab/TF_chipseq_pipeline/blob/master/README_PIPELINE.md/#troubleshooting).

### samtools ncurses bug

Prepend a directory for `libncurses.so.5` to `LD_LIBRARY_PATH`. See `install_dependencies.sh` for solution.
```
samtools: symbol lookup error: /lib/x86_64-linux-gnu/libncurses.so.5: undefined symbol: _nc_putchar
```

### Error: could not find environment: bds_atac

Unload any Anaconda Python modules. Remove locally installed Anaconda Python from your `$PATH`.

# Alternate Cloud-based Implementations

* The [Encyclopedia of DNA Elements (ENCODE) Project](https://www.encodeproject.org/pipelines/) is in the process of adopting this pipeline for uniform processing of ENCODE ATAC-seq data. The [official ENCODE implementation](https://github.com/ENCODE-DCC) by the ENCODE Data Coordination Center will be an exact mirror of our pipeline on [the DNAnexus cloud](https://www.dnanexus.com/) (i.e. results will be exactly reproducible). Note that using this service requires a user to pay for cloud compute time.

* [Epinomics](http://www.epinomics.co/) provides an independent, *free*, cloud-based pipeline implementation that adheres to the analysis protocol specifications of our pipeline. This implementation can be accessed at [https://open.epigenomics.co/#/encode](https://open.epigenomics.co/#/encode).

# Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Chuan Sheng Foo - PhD Student, Computer Science Dept., Stanford University
* Daniel Kim - MD/PhD Student, Biomedical Informatics Program, Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
