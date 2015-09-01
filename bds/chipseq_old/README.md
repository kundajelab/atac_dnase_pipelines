ENCODE ChIP-Seq Pipelines
===============================================

ENCODE ChIP-Seq pipelines are based on https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit# .

Taking advandatge of the powerful pipeline language BigDataScript (http://pcingola.github.io/BigDataScript/index.html), ENCODE ChIP-Seq pipelines have the following features:

```
1) One-command-line installation for all dependencies for ChIP-Seq pipeline
2) One command line (or one configuration file) to run the whole pipeline
3) Automatically resuming from the point of failure for failed jobs (by comparing timestamps of input/output files)
4) Optimized parallel jobs for the pipeline
5) Sun Grid Engine cluster support
6) Realtime HTML Progress report to monitor the pipeline jobs
```

### Installation instruction (pipelines and their dependencies)

```
# get the latest version of chipseq pipelines
$git clone https://github.com/kundajelab/pipelines/

# check you get the tf chipseq pipeline script
$cd pipelines/bds/chipseq
$ls -l tf_chipseq.bds

# install dependencies
$./install_dependencies.sh

# move bds.config to BigDataScript (BDS) directory
$mkdir -p $HOME/.bds
$cp bds.config $HOME/.bds/
```

Add the following lines to your $HOME/.bashrc or $HOME/.bash_profile:

```
# Java settings
export _JAVA_OPTIONS="-Xms256M -Xmx512M -XX:ParallelGCThreads=1"
export MAX_JAVA_MEM="8G"
export MALLOC_ARENA_MAX=4

# BigDataScript settings
export PATH=$PATH:$HOME/.bds
```

### Usage

There are two ways to define parameters for ChIP-Seq pipelines. For most of the parameters, they already have default values. If they are not defined in command line argument or in a configuration file, default value will be used.

1) From command line arguments 
```
# general usage
$bds tf_chipseq.bds [OPTS_FOR_PIPELINE]

# help for parameters
$bds tf_chipseq.bds -h

# example cmd. line args (human, no replicate-2 control fastq and using Anshul Kundaje's IDR)
$bds tf_chipseq.bds \
-prefix ENCSR000EGM \
-input fastq \
-fastq1 /DATA/ENCSR000EGM/ENCFF000YLW.fastq.gz \
-fastq2 /DATA/ENCSR000EGM/ENCFF000YLY.fastq.gz \
-ctl_fastq1 /DATA/ENCSR000EGM/Ctl/ENCFF000YRB.fastq.gz \
-idx_bwa /INDEX/encodeHg19Male_v0.7.10/encodeHg19Male_bwa-0.7.10.fa \
-idr_nboley false
```

2) From a configuration file
```
$bds tf_chipseq.bds [CONF_FILE]

# example configuration file (human, no replicate-2 control fastq and using Nathan Boley's IDR)
$cat [CONF_FILE]

PREFIX=ENCSR000EGM
INPUT_TYPE=fastq
INPUT_FASTQ_REP1= /DATA/ENCSR000EGM/ENCFF000YLW.fastq.gz
INPUT_FASTQ_REP2= /DATA/ENCSR000EGM/ENCFF000YLY.fastq.gz
INPUT_FASTQ_CTL_REP1= /DATA/ENCSR000EGM/Ctl/ENCFF000YRB.fastq.gz
BWA_INDEX_NAME= /INDEX/encodeHg19Male_v0.7.10/encodeHg19Male_bwa-0.7.10.fa
USE_IDR_NBOLEY=true
```

IMPORTANT! For generating bwa index, we recommend to use bwa 0.7.10.

### For cluster use (Sun Grid Engine only)

Add "-s sge" to the command line.

```
$bds -s sge tf_chipseq.bds ...
```

### Debugging and logging for BDS

```
# make BDS verbose
$bds -v tf_chipseq.bds ...

# display debugging information
$bds -d tf_chipseq.bds ...

# test run (this actually does nothing) to check input/output file names and commands
$bds -dryRun tf_chipseq.bds ...
```

For better debugging, an HTML progress report in the working directory (where you run the pipeline command) will be useful. You can monitor your BDS jobs real time.

### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
