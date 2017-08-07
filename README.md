ATAC-Seq / DNase-Seq Pipeline
===================================================

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.156534.svg)](https://doi.org/10.5281/zenodo.156534)

This pipeline is designed for automated end-to-end quality control and processing of ATAC-seq or DNase-seq data. The pipeline can be run on compute clusters with job submission engines or stand alone machines. It inherently makes uses of parallelized/distributed computing. Pipeline installation is also easy as most dependencies are automatically installed. The pipeline can be run end-to-end i.e. starting from raw FASTQ files all the way to peak calling and signal track generation; or can be started from intermediate stages as well (e.g. alignment files). The pipeline supports single-end or paired-end ATAC-seq or DNase-seq data (with or without replicates). The pipeline produces pretty HTML reports that include quality control measures specifically designed for ATAC-seq and DNase-seq data, analysis of reproducibility, stringent and relaxed thresholding of peaks, fold-enrichment and pvalue signal tracks.  The pipeline also supports detailed error reporting and easy resuming of runs. The pipeline has been tested on human, mouse and yeast ATAC-seq data and human and mouse DNase-seq data.

The ATAC-seq pipeline specification is also the official pipeline specification of the Encyclopedia of DNA Elements (ENCODE) consortium. The ATAC-seq pipeline protocol definition is [here](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit?usp=sharing). Some parts of the ATAC-seq pipeline were developed in collaboration with Jason Buenrostro, Alicia Schep and Will Greenleaf at Stanford.

The DNase-seq pipeline specification is [here](https://docs.google.com/document/d/1e3cCormg0SnQW6zr7VYBWvXC1GiPc5GSy80qlKBPwlA/edit?usp=sharing). Note that this is NOT the same as the official ENCODE DNase-seq pipeline.

* Go to [Genomic pipelines in Kundaje lab](https://kundajelab.github.io/bds_pipeline_modules)
* Go to [Discussion channel](https://groups.google.com/forum/#!forum/klab_genomic_pipelines_discuss)
* Jump to [Usage](#usage)
* Jump to [Output directory structure and file naming](#output-directory-structure-and-file-naming)
* Jump to [ENCODE accession guideline](#encode-accession-guideline)
* Jump to [Troubleshooting](#troubleshooting)

# Installation

Install software/database in a correct order according to your system. For example on Kundaje lab's clusters, you only need to install one software [Pipeline](#pipeline).

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
  * [BigDataScript](#bigdatascript)
  * [Pipeline](#pipeline)

* Stanford Sherlock cluster
  * [BigDataScript](#bigdatascript)
  * [Pipeline](#pipeline)

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

Install Miniconda3 [latest](https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh) on your system. **IMPORTANT** Make sure that the absolute path of the destination directory is short. Long path will cause an error in the depenecies installation step [issue #8](https://github.com/kundajelab/TF_chipseq_pipeline/issues/8).

```
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
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

Get the latest version of the pipeline. **Don't forget to add `--recursive`**. ATAC-Seq pipeline uses modules from an external git repo (ataqc). ATAQC will not work correctly without `--recursive`.

```
$ git clone https://github.com/kundajelab/atac_dnase_pipelines --recursive
```

## Dependencies

Install software dependencies automatically. It will create two conda environments (bds_atac and bds_atac_py3) under your conda.

```
$ bash install_dependencies.sh
```

If you don't use `install_dependencies.sh`, manually replace BDS's default `bds.config` with a correct one:

```
$ cp bds.config ./utils/bds_scr $HOME/.bds
```

If `install_dependencies.sh` fails, run `./uninstall_dependencies.sh`, fix problems and then try `bash install_dependencies.sh` again.

If you have super-user privileges on your system, it is recommended to install Miniconda3 on `/opt/miniconda3/` and share conda environment with others.

```
$ sudo su
$ bash install_dependencies.sh
$ chmod 755 -R /opt/miniconda3/  # if you get some annoying permission issues.
```

In order to make Miniconda3 accessible for all users, create an intialization script `/etc/profile.d/conda_init.sh`.

```
$ echo '#!/bin/bash' > /etc/profile.d/conda_init.sh
$ echo 'export PATH=$PATH:/opt/miniconda3/bin' >> /etc/profile.d/conda_init.sh
```

## Genome data

Install genome data for a specific genome `[GENOME]`. Currently `hg19`, `mm9`, `hg38` and `mm10` are available. Specify a directory `[DATA_DIR]` to download genome data. A species file generated on `[DATA_DIR]` will be automatically added to your `./default.env` so that the pipeline knows that you have installed genome data using `install_genome_data.sh`. If you want to install multiple genomes make sure that you use the same directory `[DATA_DIR]` for them. Each genome data will be installed on `[DATA_DIR]/[GENOME]`. If you use other BDS pipelines, it is recommended to use the same directory `[DATA_DIR]` to save disk space.

**IMPORTANT**: `install_genome_data.sh` can take longer than an hour for downloading data and building index. **DO NOT** run the script on a login node, use `qlogin` for SGE and `srun --pty bash` for SLURM.

```
# DO NOT run this on a login node
$ bash install_genome_data.sh [GENOME] [DATA_DIR]
```

If you have super-user privileges on your system, it is recommended to install genome data on `/your/data/bds_pipeline_genome_data` and share them with others.

```
$ sudo su
$ bash install_genome_data.sh [GENOME] /your/data/bds_pipeline_genome_data
```

You can find a species file `[SPECIES_FILE]` on `/your/data/bds_pipeline_genome_data` for each pipeline type. Then others can use the genome data by adding `-species_file [SPECIES_FILE_PATH]` to the pipeline command line. Or they need to add `species_file = [SPECIES_FILE_PATH]` to the section `[default]` in their `./default.env`.

# Installation for internet-free computers

The pipeline does not need internet connection but installers (`install_dependencies.sh` and `install_genome_data.sh`) do need it. So the workaround should be that first install dependencies and genome data on a computer that is connected to the internet and then move Conda and genome database directories to your internet-free one. Both computers should have **THE SAME LINUX VERSION**. 

* On your computer that has an internet access,
  * Follow [the installation instruction for general computers](#installation)
  * Move your Miniconda3 directory to `$HOME/miniconda3` on your internet-free computer.
  * Move your genome database directory, which has `bds_atac_species.conf` and directories per species, to `$HOME/genome_data` on your internet-free computer. `$HOME/genome_data` on your internet-free computer should have `bds_atac_species.conf`.
  * Move your BDS directory `$HOME/.bds` to `$HOME/.bds` on your internet-free computer.
  * Move your pipeline directory `atac_dnase_pipelines/` to `$HOME/atac_dnase_pipelines/` on your internet-free computer.

* On your internet-free computer,
  * Add your `miniconda3/bin` and BDS binary to `$PATH` in your bash initialization script (`$HOME/.bashrc` or `$HOME/.bash_profile`).

     ```
     export PATH="$PATH:$HOME/miniconda3/bin"
     export PATH="$PATH:$HOME/.bds"
     ```

  * Modify `[default]` section in `$HOME/atac_dnase_pipelines/default.env`.

     ```
     [default]
     conda_bin_dir=$HOME/miniconda3/bin
     species_file=$HOME/genome_data/bds_atac_species.conf
     ```

* Modify all paths in `$HOME/genome_data/bds_atac_species.conf` so that they correctly point to the right files.
* Check BDS version.
     ```
     $ bds -version
     Bds 0.99999e (build 2016-08-26 06:34), by Pablo Cingolani
     ```
* Make sure that your java rumtime version is >= 1.8.
     ```
     $ java -version
     java version "1.8.0_111"
     Java(TM) SE Runtime Environment (build 1.8.0_111-b14)
     Java HotSpot(TM) 64-Bit Server VM (build 25.111-b14, mixed mode)
     ```

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
IDR analysis is optional in the pipeline by default. If there are more than two replicates, IDR will be done for every pair of replicates. to enable IDR add the following:
```
-enable_idr
```
For multimapping, 
```
-multimapping [NO_MULTIMAPPING; 4 by default]
```
To force a set of parameters (`-smooth_win 73 -idr_thresh 0.05 -multimapping 4`) for ENCODE3.
```
-ENCODE3
```

You can also define parameters in a configuration file. Key names in a configruation file are identical to parameter names on command line. 
```
$ bds atac.bds [CONF_FILE]

$ cat [CONF_FILE]
species = [SPECIES; hg19, mm9, ...]
nth   = [NUM_THREADS]
fastq1_1= [READ1]
fastq1_2= [READ2]
...
```

## List of all parameters

To list all parameters and default values for them,
```
$ bds atac.bds

== atac pipeline settings
        -type <string>                   : Type of the pipeline. atac-seq or dnase-seq (default: atac-seq).
        -dnase_seq <bool>                : DNase-Seq (no tn5 shifting).
        -trimmed_fastq <bool>            : Skip fastq-trimming stage.
        -align <bool>                    : Align only (no MACS2 peak calling or IDR or ataqc analysis).
        -subsample_xcor <string>         : # reads to subsample for cross corr. analysis (default: 25M).
        -subsample <string>              : # reads to subsample exp. replicates. Subsampled tagalign will be used for steps downstream (default: 0; no subsampling).
        -true_rep <bool>                 : No pseudo-replicates.
        -no_ataqc <bool>                 : No ATAQC
        -no_xcor <bool>                  : No Cross-correlation analysis.
        -csem <bool>                     : Use CSEM for alignment.
        -smooth_win <string>             : Smoothing window size for MACS2 peak calling (default: 150).
        -idr_thresh <real>               : IDR threshold : -log_10(score) (default: 0.1).
        -old_trimmer <bool>              : Use legacy trim adapters (trim_galore and trimAdapter.py).
        -ENCODE3 <bool>                  : Force to use parameter set (-smooth_win 73 -idr_thresh 0.05 -multimapping 4) for ENCODE3.
        -ENCODE <bool>                   : Force to use parameter set (-smooth_win 73 -idr_thresh 0.05 -multimapping 4) for ENCODE.
        -no_browser_tracks <bool>        : Disable generation of genome browser tracks (workaround for bzip2 shared library issue).
        -overlap_pval_thresh <real>      : p-val threshold for overlapped peaks (default: 0.01).
        -macs2_pval_thresh <real>        : MACS2 p-val threshold for calling peaks (default: 0.1).
        -macs2_pval_thresh_bw <real>     : MACS2 p-val threshold for generating BIGWIG signal tracks (default: 0.1).
        -enable_idr <bool>               : Enable IDR on called peaks.
== configuration file settings
        -c <string>                      : Configuration file path.
        -env <string>                    : Environment file path.
== parallelization settings
        -no_par <bool>                   : Serialize all tasks (individual tasks can still use multiple threads up to '-nth').
        -nth <int>                       : Maximum # threads for a pipeline. (default: 8).
== cluster/system/resource settings
        -wt <string>                     : Walltime for all single-threaded tasks (example: 8:10:00, 3h, 3600, default: 5h50m, 5:50:00).
        -memory <string>                 : Maximum memory for all single-threaded tasks (equivalent to '-mem', example: 4.5G, 1024M, default: 7G).
        -use_system <string>             : Force to use a system (equivalent to 'bds -s [SYSTEM_NAME] ...', any system defined in bds.config can be used).
        -nice <int>                      : Set process priority for all tasks (default: 0; -20 (highest) ~ 19 (lowest) ).
        -retrial <int>                   : # of Retrial for failed tasks (default: 0).
        -q <string>                      : Submit tasks to a specified cluster queue.
        -unlimited_mem_wt <bool>         : Use unlimited max. memory and walltime.
== shell environment settings
        -mod <string>                    : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
        -shcmd <string>                  : Shell commands separated by ;. Shell var. must be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test").
        -addpath <string>                : Path separated by ; or : to be PREPENDED to \$PATH (example: "/bin/test:${HOME}/utils").
        -conda_env <string>              : Anaconda Python (or Miniconda) environment name for all softwares including Python2.
        -conda_env_py3 <string>          : Anaconda Python (or Miniconda) environment name for Python3.
        -conda_bin_dir <string>          : Anaconda Python (or Miniconda) bin directory.
        -cluster_task_min_len <int>      : Minimum length for a cluster job in seconds (dealing with NFS delayed write, default: 60).
        -cluster_task_delay <int>        : Constant delay for every job in seconds (dealing with NFS delayed write, default: 0).
== output/title settings
        -out_dir <string>                : Output directory (default: out).
        -title <string>                  : Prefix for HTML report and outputs without given prefix.
== species settings
        -species <string>                : Species. need to specify '-species_file' too if you have not installed genome database with 'install_genome_data.sh'.
        -species_file <string>           : Species file path.
        -species_browser <string>        : Species name in WashU genome browser.
        -ref_fa <string>                 : Reference genome sequence fasta.
        -chrsz <string>                  : Chromosome sizes file path (use fetchChromSizes from UCSC tools).
        -blacklist <string>              : Blacklist bed.
        -seq_dir <string>                : Reference genome sequence directory path (where chr*.fa exist).
== ENCODE accession settings
        -ENCODE_accession <string>       : ENCODE experiment accession ID (or dataset).
        -ENCODE_award_rfa <string>       : ENCODE award RFA (e.g. ENCODE3).
        -ENCODE_assay_category <string>  : ENCODE assay category.
        -ENCODE_assay_title <string>     : ENCODE assay title.
        -ENCODE_award <string>           : ENCODE award (e.g. /awards/U41HG007000/).
        -ENCODE_lab <string>             : Lab (e.g. /labs/anshul-kundaje/)
        -ENCODE_assembly <string>        : hg19, GRCh38, mm9, mm10.
        -ENCODE_alias_prefix <string>    : Alias prefix, Alias = alias_prefix: + filename + alias_suffix
        -ENCODE_alias_suffix <string>    : Alias suffix, Alias = alias_prefix: + filename + alias_suffix
== report settings
        -url_base <string>               : URL base for output directory.
        -viz_genome_coord <string>       : WashU genome browser genome coordinate (e.g. chr7:27117661-27153380).
== fastq input definition :
        Single-ended : For replicate '-fastq[REP_ID]', For control '-ctl_fastq[REP_ID]'
        Paired end : For replicate '-fastq[REP_ID]_[PAIR_ID]', For control '-ctl_fastq[REP_ID]_[PAIR_ID]'
== bam input (raw or filtered) definition :
        Raw bam : For replicate '-bam[REP_ID]', For control '-ctl_bam[REP_ID]'.
        Filtered bam : For replicate '-filt_bam[REP_ID]', For control '-ctl_filt_bam[REP_ID]'.
== tagalign input definition :
        For replicate '-tag[REP_ID]', For control '-ctl_tag[REP_ID]'.
== narrow peak input definition :
        For true replicates, use '-peak1' and '-peak2',
        For pooled replicates, use '-peak_pooled',
        For two PR (self-pseudo-replicates), use '-peak[REP_ID]_pr1' and '-peak[REP_ID]_pr2'
        For two PPR (pooled pseudo-replicates), use '-peak_ppr1' and '-peak_ppr2'
== input endedness settings (SE or PE) :
        -se <bool>                       : Singled-ended data set. To specify it for each replicate, '-se[REP_ID]' for exp. reps, '-ctl_se[CTL_ID]' for control.
        -pe <bool>                       : Paired end data set. To specify it for each replicate, '-pe[REP_ID]' for exp. reps, '-ctl_pe[CTL_ID]' for controls.
== adapter sequence definition :
        Single-ended : For replicate '-adapter[REP_ID]'
        Paired end : For replicate '-adapter[REP_ID]_[PAIR_ID]'
== align multimapping settings
        -multimapping <int>              : # alignments reported for multimapping (default: 0).
== align bowtie2 settings (requirements: -bwt2_idx)
        -bwt2_idx <string>               : Bowtie2 index (full path prefix of *.1.bt2 file).
        -scoremin_bwt2 <string>          : Replacement --score-min for bowtie2.
        -wt_bwt2 <string>                : Walltime for bowtie2 (default: 47h, 47:00:00).
        -mem_bwt2 <string>               : Max. memory for bowtie2 (default: 12G).
        -extra_param_bwt2 <string>       : Extra parameter for bowtie2.
== adapter trimmer settings
        -adapter_err_rate <string>       : Maximum allowed adapter error rate (# errors divided by the length of the matching adapter region, default: 0.10).
        -min_trim_len <int>              : Minimum trim length for cutadapt -m, throwing away processed reads shorter than this (default: 5).
        -wt_trim <string>                : Walltime for adapter trimming (default: 23h, 23:00:00).
        -mem_trim <string>               : Max. memory for adapter trimming (default: 12G).
== postalign bam settings
        -mapq_thresh <int>               : Threshold for low MAPQ reads removal (default: 30).
        -rm_chr_from_tag <string>        : Perl style reg-ex to exclude reads from tag-aligns. (example: 'other|ribo|mito|_', '_', default: blank)
        -no_dup_removal <bool>           : No dupe removal when filtering raw bam.
        -wt_dedup <string>               : Walltime for post-alignment filtering (default: 23h, 24:00:00).
        -mem_dedup <string>              : Max. memory for post-alignment filtering (default: 12G).
        -dup_marker <string>             : Dup marker for filtering mapped reaads in BAMs: picard or sambamba (default: picard).
        -use_sambamba_markdup <bool>     : Use sambamba markdup instead of Picard MarkDuplicates (default: false).
== postalign bed/tagalign settings
        -mem_shuf <string>               : Max. memory for UNIX shuf (default: 12G).
        -fraglen0 <bool>                 : (LEGACY PARAM) Set predefined fragment length as zero for cross corr. analysis (add -speak=0 to run_spp.R).
        -speak_xcor <int>                : Set user-defined cross-corr. peak strandshift (-speak= in run_spp.R). Use -1 to disable (default: -1).
        -extra_param_xcor <string>       : Set extra parameters for run_spp.R (cross-corr. analysis only).
== postalign bed/tagalign settings
        -mem_shuf <string>               : Max. memory for UNIX shuf (default: 12G).
        -fraglen0 <bool>                 : (LEGACY PARAM) Set predefined fragment length as zero for cross corr. analysis (add -speak=0 to run_spp.R).
        -speak_xcor <int>                : Set user-defined cross-corr. peak strandshift (-speak= in run_spp.R). Use -1 to disable (default: -1).
        -extra_param_xcor <string>       : Set extra parameters for run_spp.R (cross-corr. analysis only).
== callpeak macs2 settings (requirements: -chrsz -gensz)
        -gensz <string>                  : Genome size; hs for human, mm for mouse.
        -wt_macs2 <string>               : Walltime for MACS2 (default: 23h, 23:00:00).
        -mem_macs2 <string>              : Max. memory for MACS2 (default: 15G).
        -extra_param_macs2 <string>      : Extra parameters for macs2 callpeak.
== callpeak naive overlap settings
        -nonamecheck <bool>              : bedtools intersect -nonamecheck (bedtools>=2.24.0, use this if you get bedtools intersect naming convenction warnings/errors).
== callpeak etc settings
        -npeak_filt <int>                : # top peaks filtered from a narrow peak files (default: 500000).
== IDR settings
        -idr_suffix <bool>               : Append IDR threshold to IDR output directory.
== ATAQC settings
        -tss_enrich <string>             : TSS enrichment bed for ataqc.
        -dnase <string>                  : DNase bed (open chromatin region file) for ataqc.
        -prom <string>                   : Promoter bed (promoter region file) for ataqc.
        -enh <string>                    : Enhancer bed (enhancer region file) for ataqc.
        -reg2map <string>                : Reg2map (file with cell type signals) for ataqc.
        -reg2map_bed <string>            : Reg2map_bed (file of regions used to generate reg2map signals) for ataqc.
        -roadmap_meta <string>           : Roadmap metadata for ataqc.
        -mem_ataqc <string>              : Max. memory for ATAQC (default: 20G).
        -wt_ataqc <string>               : Walltime for ATAQC (default: 47h, 47:00:00).
```

## Stopping / Resuming pipeline

Press Ctrl + C on a terminal or send any kind of kill signals to it. Please note that this will delete all intermediate files and incomplete outputs for the running tasks. The pipeline automatically determines if each task has finished or not (by comparing timestamps of input/output files for each task). To run the pipeline from the point of failure, correct error first and then just run the pipeline with the same command that you started the pipeline with. There is no additional parameter for restarting the pipeline.

## Running pipelines with a cluster engine

On servers with a cluster engine (such as Sun Grid Engine and SLURM), **DO NOT QSUB/SBATCH BDS COMMAND LINE**. Run BDS command directly on login nodes. BDS is a task manager and it will automatically submit(qsub/sbatch) and manage its sub tasks.

**IMPORTANT!** Please read this section carefully if you run pipelines on Stanford SCG4 and Sherlock cluster.

Most clusters have a policy to limit number of threads and memory per user on a login node. One BDS process, as a Java-based task manager, takes up to 1GB of memory and 50 threads even though it just submits/monitors subtasks. So if you want to run more than 50 pipelines in parallel, your cluster will kill BDS processes due to resource limit on a login node (check resource limit per user with `ulimit -a`). For example of 50 pipelines, 50 GB of memory and 2500 threads will be taken by 50 BDS processes. So the Workaround for this is to make an interactive node to keep all BDS processes alive. Such interactive node must have long walltime enough to wait for all pipelines in it to finish. Recommended resource setting is 0.5GB memory per pipeline.

SGE example to make an interactive node for 100 pipelines: 1 cpu, 100GB memory, 3 days walltime.

```
$ qlogin -l h_rt=72:00:00 -l h_vmem=100G
```

SLURM example to make an interactive node for 100 pipelines: 1 cpus, 100GB memory, 3 days walltime.

```
$ srun -n 1 --mem 100G -t 3-0 -p [YOUR_PARTITON] --qos normal --pty /bin/bash -i -l 
```

Once you get an interactive node, repeat the following commands per sample to run a pipeline with using [`bds_scr`](#managing-multiple-pipelines).

```
$ cd [WORK_DIR]
$ bds_scr [SCREEN_NAME] [LOG_FILE_PATH] atac.bds -q [SGE_QUEUE_OR_SLURM_PARTITION] -nth [MAX_NUM_THREAD_PER_PIPELINE] ...
$ sleep 2 # wait for 2 seconds for safety
```

Then you can monitor your pipelines with `screen -ls` and `tail -f [LOG_FILE_PATH]`. If you want to run more than 200 pipelines, you would want to make multiple interactive nodes and distribute your samples to them.

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
$ bds_scr [SCR_NAME] [LOG_FILE_NAME] atac.bds ...
$ bds_scr [SCR_NAME] atac.bds ...
$ bds_scr [SCR_NAME] -s sge atac.bds ...
```

Once the pipeline run is done, the screen will be automatically closed. To kill a pipeline manually while it's running, use `./utils/kill_scr` or `screen -X quit`:

```
$ screen -X -S [SCR_NAME] quit
$ kill_scr [SCR_NAME]
```

## Java issues (memory and temporary directory)

Picard tools is used for marking dupes in the reads and it's based on Java. If you see any Java heap space errors then increase memory limit for Java with `-mem_ataqc [MEM]` (default: `20G`) and `-mem_dedup [MEM]` (default: `12G`).

If your `/tmp` quickly fills up and you want to change temporary directory for all Java apps in the pipeline, then add the following line to your bash startup script (`$HOME/.bashrc`). Our pipeline takes in `$TMPDIR` (not `$TMP`) for all Java apps.

```
export TMPDIR=/your/temp/dir/
```

Another quick workaround for dealing with Java issues is not to use Picard tools in the pipeline. Add `-use_sambamba_markdup` to your command line and then you can use `sambamba markdup` instead of `picard markdup`.


## How to customize genome data installer?

Please refer to the section `Installer for genome data` on [BDS pipeline programming](https://kundajelab.github.io/bds_pipeline_modules/programming.html).


## Useful HTML reports

There are two kinds of HTML reports provided by the pipeline.

* BigDataScript HTML report for debugging: Located at the working folder with name atac_[TIMESTAMP]_report.html. This report is automatically generated by BigDataScript and is useful for debugging since it shows summary, timeline, Stdout and Stderr for each job.

* ATAC-Seq pipeline report for QC and result: The pipeline automatically generate a nice HTML report (Report.html) on its output directory (specified with -out_dir or just './out'). It summarizes files and directory structure, includes QC reports and show a workflow diagram and genome browser tracks for peaks and signals (bigwigs for pValue and fold change). Move your output directory to a web directory (for example, /var/www/somewhere) or make a softlink of it to a web directory. For genome browser tracks, specify your web directory root for your output  While keeping its structure. Make sure that you have bgzip and tabix installed on your system. Add the following to the command line:

      -url_base http://your/url/to/output -title [PREFIX_FOR_YOUR_REPORT]


# Output directory structure and file naming

For more details, refer to the file table section in an HTML report generated by the pipeline. Files marked as (E) are outputs to be uploaded during ENCODE accession.
```
out                               # root dir. of outputs
│
├ *report.html                    #  HTML report
├ *tracks.json                    #  Tracks datahub (JSON) for WashU browser
├ ENCODE_summary.json             #  Metadata of all datafiles and QC results
│
├ align                           #  mapped alignments
│ ├ rep1                          #   for true replicate 1 
│ │ ├ *.trim.fastq.gz             #    adapter-trimmed fastq
│ │ ├ *.bam                       #    raw bam
│ │ ├ *.nodup.bam (E)             #    filtered and deduped bam
│ │ ├ *.tagAlign.gz               #    tagAlign (bed6) generated from filtered bam
│ │ ├ *.tn5.tagAlign.gz           #    TN5 shifted tagAlign for ATAC pipeline (not for DNase pipeline)
│ │ └ *.*M.tagAlign.gz            #    subsampled tagAlign for cross-corr. analysis
│ ├ rep2                          #   for true repilicate 2
│ ...
│ ├ pooled_rep                    #   for pooled replicate
│ ├ pseudo_reps                   #   for self pseudo replicates
│ │ ├ rep1                        #    for replicate 1
│ │ │ ├ pr1                       #     for self pseudo replicate 1 of replicate 1
│ │ │ ├ pr2                       #     for self pseudo replicate 2 of replicate 1
│ │ ├ rep2                        #    for repilicate 2
│ │ ...                           
│ └ pooled_pseudo_reps            #   for pooled pseudo replicates
│   ├ ppr1                        #    for pooled pseudo replicate 1 (rep1-pr1 + rep2-pr1 + ...)
│   └ ppr2                        #    for pooled pseudo replicate 2 (rep1-pr2 + rep2-pr2 + ...)
│
├ peak                             #  peaks called
│ └ macs2                          #   peaks generated by MACS2
│   ├ rep1                         #    for replicate 1
│   │ ├ *.narrowPeak.gz            #     narrowPeak (p-val threshold = 0.01)
│   │ ├ *.filt.narrowPeak.gz (E)   #     blacklist filtered narrowPeak 
│   │ ├ *.narrowPeak.bb (E)        #     narrowPeak bigBed
│   │ ├ *.narrowPeak.hammock.gz    #     narrowPeak track for WashU browser
│   │ ├ *.pval0.1.narrowPeak.gz    #     narrowPeak (p-val threshold = 0.1)
│   │ └ *.pval0.1.*K.narrowPeak.gz #     narrowPeak (p-val threshold = 0.1) with top *K peaks
│   ├ rep2                         #    for replicate 2
│   ...
│   ├ pseudo_reps                          #   for self pseudo replicates
│   ├ pooled_pseudo_reps                   #   for pooled pseudo replicates
│   ├ overlap                              #   naive-overlapped peaks
│   │ ├ *.naive_overlap.narrowPeak.gz      #     naive-overlapped peak
│   │ └ *.naive_overlap.filt.narrowPeak.gz #     naive-overlapped peak after blacklist filtering
│   └ idr                           #   IDR thresholded peaks
│     ├ true_reps                   #    for replicate 1
│     │ ├ *.narrowPeak.gz           #     IDR thresholded narrowPeak
│     │ ├ *.filt.narrowPeak.gz (E)  #     IDR thresholded narrowPeak (blacklist filtered)
│     │ └ *.12-col.bed.gz           #     IDR thresholded narrowPeak track for WashU browser
│     ├ pseudo_reps                 #    for self pseudo replicates
│     │ ├ rep1                      #    for replicate 1
│     │ ...
│     ├ optimal_set                 #    optimal IDR thresholded peaks
│     │ └ *.filt.narrowPeak.gz (E)  #     IDR thresholded narrowPeak (blacklist filtered)
│     ├ conservative_set            #    optimal IDR thresholded peaks
│     │ └ *.filt.narrowPeak.gz (E)  #     IDR thresholded narrowPeak (blacklist filtered)
│     ├ pseudo_reps                 #    for self pseudo replicates
│     └ pooled_pseudo_reps          #    for pooled pseudo replicate
│
│   
│ 
├ qc                              #  QC logs
│ ├ *IDR_final.qc                 #   Final IDR QC
│ ├ rep1                          #   for true replicate 1
│ │ ├ *.align.log                 #    Bowtie2 mapping stat log
│ │ ├ *.dup.qc                    #    Picard (or sambamba) MarkDuplicate QC log
│ │ ├ *.pbc.qc                    #    PBC QC
│ │ ├ *.nodup.flagstat.qc         #    Flagstat QC for filtered bam
│ │ ├ *M.cc.qc                    #    Cross-correlation analysis score for tagAlign
│ │ ├ *M.cc.plot.pdf/png          #    Cross-correlation analysis plot for tagAlign
│ │ └ *_qc.html/txt               #    ATAQC report
│ ...
│
├ signal                          #  signal tracks
│ ├ macs2                         #   signal tracks generated by MACS2
│ │ ├ rep1                        #    for true replicate 1 
│ │ │ ├ *.pval.signal.bigwig (E)  #     signal track for p-val
│ │ │ └ *.fc.signal.bigwig   (E)  #     signal track for fold change
│ ...
│ └ pooled_rep                    #   for pooled replicate
│ 
└ report                          # files for HTML report
```

# ENCODE accession guideline

For each pipeline rune, `ENCODE_summary.json` file is generated under the output directory (`-out_dir`) for ENCODE accession (uploading pipeline outputs to the ENCODE portal). This JSON file includes all metadata and QC metrics required for ENCODE accession.

For ENCODE3, Please make sure that you run pipelines with `-ENCODE3` flag.

Parameters required for ENCODE accesssion:
```
# required
        -ENCODE_accession <string>       : ENCODE experiment accession ID (or dataset).
        -ENCODE_award <string>           : ENCODE award (e.g. /awards/U41HG007000/).
        -ENCODE_lab <string>             : Lab (e.g. /labs/anshul-kundaje/)
        -ENCODE_assembly <string>        : hg19, GRCh38, mm9, mm10.
        -ENCODE_alias_prefix <string>    : Alias = Alias_prefix + filename
# optional
        -ENCODE_award_rfa <string>       : ENCODE award RFA (e.g. ENCODE3).
        -ENCODE_assay_category <string>  : ENCODE assay category.
        -ENCODE_assay_title <string>     : ENCODE assay title.
```

We also provide an [ENCODE fastq downloader](https://github.com/kundajelab/ENCODE_downloader). It downloads fastqs matching award_rfa, assay_category and assay_title, and then automatically generate a shell script to run multiple pipelines. Such shell script also includes these ENCODE accession parameter set.

## ENCODE accession spreadsheet (CSV) generation

`./utils/parse_summary_ENCODE_accession_recursively.py` recursively finds `ENCODE_summary.json` files and parse them to generate one big CSV spreadsheet for ENCODE accession.

```
$ python ./utils/parse_summary_ENCODE_accession_recursively.py -h

usage: ENCODE_summary.json parser for ENCODE accession [-h]
                                                       [--out-file OUT_FILE]
                                                       [--search-dir SEARCH_DIR]
                                                       [--json-file JSON_FILE]
                                                       [--sort-by-genome-and-exp]
                                                       [--ignored-accession-ids-file IGNORED_ACCESSION_IDS_FILE]

Recursively find ENCODE_summary.json, parse it and make a CSV for uploading to
the ENCODE portal. Use https://github.com/ENCODE-DCC/pyencoded-
tools/blob/master/ENCODE_submit_files.py for uploading.

optional arguments:
  -h, --help            show this help message and exit
  --out-file OUT_FILE   Output CSV filename)
  --search-dir SEARCH_DIR
                        Root directory to search for ENCODE_summary.json
  --json-file JSON_FILE
                        Specify json file name to be parsed
  --sort-by-genome-and-exp
                        Sort rows by genomes and ENCODE experiment accession
                        ID
  --ignored-accession-ids-file IGNORED_ACCESSION_IDS_FILE
                        Accession IDs in this text file will be ignored. (1
                        acc. ID per line)
```

## QC metrics spreadsheet (TSV) generation

`./utils/parse_summary_qc_recursively.py` recursively finds `ENCODE_summary.json` files and parse them to generate one big TSV spreadsheet for QC metrics.

```
$ python parse_summary_qc_recursively.py -h
usage: ENCODE_summary.json parser for QC [-h] [--out-file OUT_FILE]
                                         [--search-dir SEARCH_DIR]
                                         [--json-file JSON_FILE]

Recursively find ENCODE_summary.json, parse it and make a TSV spreadsheet of
QC metrics.

optional arguments:
  -h, --help            show this help message and exit
  --out-file OUT_FILE   Output TSV filename)
  --search-dir SEARCH_DIR
                        Root directory to search for ENCODE_summary.json
  --json-file JSON_FILE
                        Specify json file name to be parsed
```


# Programming with BDS

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

### Error: could not find environment: bds_atac

Unload any Anaconda Python modules. Remove locally installed Anaconda Python from your `$PATH`.

# Alternate Cloud-based Implementations

* The [Encyclopedia of DNA Elements (ENCODE) Project](https://www.encodeproject.org/pipelines/) is in the process of adopting this pipeline for uniform processing of ENCODE ATAC-seq data. The [official ENCODE implementation](https://github.com/ENCODE-DCC) by the ENCODE Data Coordination Center will be an exact mirror of our pipeline on [the DNAnexus cloud](https://www.dnanexus.com/) (i.e. results will be exactly reproducible). Note that using this service requires a user to pay for cloud compute time.

* [Epinomics](http://www.epinomics.co/) provides an independent, *free*, cloud-based pipeline implementation that adheres to the analysis protocol specifications of our pipeline. This implementation can be accessed at [https://open.epigenomics.co/#/encode](https://open.epigenomics.co/#/encode).

# Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Chuan Sheng Foo - PhD Student, Computer Science Dept., Stanford University
* Daniel Kim - MD/PhD Student, Biomedical Informatics Program, Stanford University
* Nathan Boley - Postdoc, Dept. of Genetics, Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University

We'd also like to acknowledge Jason Buenrostro, Alicia Schep and William Greenleaf who contributed prototype code for some parts of the ATAC-seq pipeline.
