ATAC-Seq / DNase-Seq Pipeline
===================================================

* Go to [Discussion channel](https://groups.google.com/forum/#!forum/klab_genomic_pipelines_discuss)
* Jump to [Usage](#usage)
* Jump to [Troubleshooting](#troubleshooting)

# Installation

* General computer
  * [Java](#java)
  * [Conda](#conda)
  * [BigDataScript](#bigdatascript)
  * [Pipeline](#pipeline)
  * [Dependencies](#dependencies)
  * [Genome data](#genome-data)

* Kundaje lab's clusters
  * [Pipeline](#pipeline)

* Stanford SCG cluster
  * [Conda](#conda)
  * [BigDataScript](#bigdatascript)
  * [Pipeline](#pipeline)
  * [Dependencies](#dependencies)

* Stanford Sherlock cluster
  * [Conda](#conda)
  * [BigDataScript](#bigdatascript)
  * [Pipeline](#pipeline)
  * [Dependencies](#dependencies)

## Java

Install Java 8 (jdk >= 1.8 or jre >= 1.8) on your system. If you don't have super-user privileges on your system, locally install it and add it to your `$PATH`.

* For Debian/Ubuntu (>14.10) based Linux:

     ```
     $ sudo apt-get install git openjdk-8-jre
     ```

* For Fedora/Red-Hat based Linux: 
 
     ```
     $ sudo yum install git java-1.8.0-openjdk
     ```

* For Ubuntu 14.04 (trusty):

     ```
     $ sudo add-apt-repository ppa:webupd8team/java -y
     $ sudo apt-get update
     $ sudo apt-get install oracle-java8-installer
     ```

## Conda

Install Miniconda3 [4.0.5](https://repo.continuum.io/miniconda/Miniconda3-4.0.5-Linux-x86_64.sh) on your system. Recent versions of conda (>4.0.10) is buggy in parallel activation and do not work correctly with the pipeline. If you already have your own conda, downgrade it to 4.0.5 (`conda install conda==4.0.5`).

```
$ wget https://repo.continuum.io/miniconda/Miniconda3-4.0.5-Linux-x86_64.sh
$ bash Miniconda3-4.0.5-Linux-x86_64.sh
```

Answer `yes` for the final question. If you choose `no`, you need to manually add Miniconda3 to your `$HOME/.bashrc`.

```
Do you wish the installer to prepend the Miniconda3 install location
to PATH in your /your/home/.bashrc ? [yes|no]
[no] >>> yes
```

Remove any other Anaconda from your `$PATH`. Check your loaded modules with `$ module list` and unload any Anaconda modules in your `$HOME/.bashrc`. Open a new terminal after installation.

## BigDataScript

Install BigDataScript v0.99999e ([forked](https://github.com/leepc12/BigDataScript)) on your system.
The original [BDS v0.99999e](https://github.com/pcingola/BigDataScript) does not work correctly with the pipeline
(see [PR #142](https://github.com/pcingola/BigDataScript/pull/142) and [issue #131](https://github.com/pcingola/BigDataScript/issues/131)).

```
$ wget https://github.com/leepc12/BigDataScript/blob/master/distro/bds_Linux.tgz?raw=true -O bds_Linux.tgz
$ mv bds_Linux.tgz $HOME && cd $HOME && tar zxvf bds_Linux.tgz
```

Add `export PATH=$PATH:$HOME/.bds` to your `$HOME/.bashrc`. If Java memory occurs, add `export _JAVA_OPTIONS="-Xms256M -Xmx728M -XX:ParallelGCThreads=1"` too.


## Pipeline

Get the latest version of the pipeline. **Don't forget to add `--recursive`**. ATAC-Seq pipeline uses modules in two external git repos (ataqc, TF_chipseq_pipeline). It will not work correctly without `--recursive`.

```
$ git clone https://github.com/kundajelab/bds_atac --recursive
```

## Dependencies

Install software dependencies automatically. It will create two conda environments (bds_atac and bds_atac_py3) under your conda.

```
$ ./install_dependencies.sh
```

If you see the following error, see [issue #8](https://github.com/kundajelab/TF_chipseq_pipeline/issues/8)

```
Error: ERROR: placeholder '/root/miniconda3/envs/_build_placehold_placehold_placehold_placehold_placehold_p' too short in: glib-2.43.0-2
```

If you don't use `install_dependencies.sh`, manually replace BDS's default `bds.config` with a correct one:

```
$ cp bds.config ./utils/bds_scr $HOME/.bds
```

If `install_dependencies.sh` fails, run `./uninstall_dependencies.sh`, fix problems and then try `./install_dependencies.sh` again.

If you have super-user privileges on your system, it is recommended to install Miniconda3 on `/opt/miniconda3/` and share conda environment with others.

```
$ sudo su
$ ./install_dependencies.sh
$ chmod 755 -R /opt/miniconda3/  # if you get some annoying permission issues.
```

In order to make Miniconda3 accessible for all users, create an intialization script `/etc/profile.d/conda_init.sh`.

```
$ echo '#!/bin/bash' > /etc/profile.d/conda_init.sh
$ echo 'export PATH=$PATH:/opt/miniconda3/bin' >> /etc/profile.d/conda_init.sh
```

## Genome data

Install genome data for a specific genome `[GENOME]`. Currently `hg19`, `mm9`, `hg38`(BETA) and `mm10`(BETA) are available. Specify a directory `[DATA_DIR]` to download genome data. A species file generated on `[DATA_DIR]` will be automatically added to your `./default.env` so that the pipeline knows that you have installed genome data using `install_genome_data.sh`. If you want to install multiple genomes make sure that you use the same directory `[DATA_DIR]` for them. Each genome data will be installed on `[DATA_DIR]/[GENOME]`. If you use other BDS pipelines, it is recommended to use the same directory `[DATA_DIR]` to save disk space.

**IMPORTANT**: `install_genome_data.sh` can take longer than an hour for downloading data and building index. **DO NOT** run the script on a login node, use `qlogin` for SGE and `sdev` for SLURM.

```
# DO NOT run this on a login node
$ ./install_genome_data.sh [GENOME] [DATA_DIR]
```

If you have super-user privileges on your system, it is recommended to install genome data on `/your/data/bds_pipeline_genome_data` and share them with others.

```
$ sudo su
$ ./install_genome_data.sh [GENOME] /your/data/bds_pipeline_genome_data
```

You can find a species file `[SPECIES_FILE]` on `/your/data/bds_pipeline_genome_data` for each pipeline type. Then others can use the genome data by adding `-species_file [SPECIES_FILE_PATH]` to the pipeline command line. Or they need to add `species_file = [SPECIES_FILE_PATH]` to the section `[default]` in their `./default.env`.

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

The pipeline can automatically detect adapters if you remove `-adapter` from your command line.

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
```
-se1 -pe2 -se3 ...
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

## Specifying cluster queue

You can let BDS submit its subtasks to a specific queue `[QUEUE_NAME]` on Sun Grid Engine.
```
bds -q [QUEUE_NAME] -s sge atac.bds ...
bds -s sge atac.bds -q [QUEUE_NAME] ...
```

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

## Programming with BDS

* [Using genomic pipeline modules in Kundaje lab](https://kundajelab.github.io/bds_pipeline_modules/programming.html)

* [BigDataScript github repo](https://github.com/pcingola/BigDataScript)

* [BigDataScript documentation](http://pcingola.github.io/BigDataScript/bigDataScript_manual.html)

# Requirements 

* For python2 (python 2.x >= 2.7) and R-3.x, [requirements.txt](requirements.txt)

* For python3, [requirements_py3.txt](requirements_py3.txt)

# Troubleshooting

See more [troubleshootings](https://kundajelab.github.io/bds_pipeline_modules/troubleshooting.html)

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
