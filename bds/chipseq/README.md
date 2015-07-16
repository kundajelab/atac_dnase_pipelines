ChIP-Seq Pipelines
===============================================

ChIP-Seq pipelines are based on https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#.

### Installation instruction (pipelines and their dependencies)

```
# get the latest version of chipseq pipelines
$git clone https://github.com/kundajelab/pipelines/

# tf chipseq pipeline script
$ls -l pipelines/bds/chipseq/tf_chipseq.bds

# move bds.config to BigDataScript (BDS) directory
$mkdir -p $HOME/.bds
$cp pipelines/bds/bds.config $HOME/.bds/

# install dependencies
$cd pipelines/bds/chipseq
$./install_dependencies.sh
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

# example cmd. line args (no replicate-2 control fastq and using Anshul Kundaje's IDR)
$bds tf_chipseq.bds \
-prefix ENCSR000EGM \
-input fastq \
-fastq1 /DATA/ENCSR000EGM/ENCFF000YLW.fastq.gz \
-fastq2 /DATA/ENCSR000EGM/ENCFF000YLY.fastq.gz \
-ctl_fastq1_1 /DATA/ENCSR000EGM/Ctl/ENCFF000YRB.fastq.gz \
-idx_bwa /INDEX/encodeHg19Male_v0.7.10/encodeHg19Male_bwa-0.7.10.fa \
-idr_nboley false
```

2) From a configuration file
```
$bds tf_chipseq.bds [CONF_FILE]

# example configuration file (no replicate-2 control fastq and using Nathan Boley's IDR)
$cat [CONF_FILE]

PREFIX=ENCSR000EGM
INPUT_TYPE=fastq
INPUT_FASTQ_REP1= /srv/scratch/leepc12/data/ENCODE/ENCSR000EGM/ENCFF000YLW.fastq.gz
INPUT_FASTQ_REP2= /srv/scratch/leepc12/data/ENCODE/ENCSR000EGM/ENCFF000YLY.fastq.gz
INPUT_FASTQ_CTL_REP1= /srv/scratch/leepc12/data/ENCODE/ENCSR000EGM/Ctl/ENCFF000YRB.fastq.gz
BWA_INDEX_NAME= /srv/scratch/leepc12/hg19/encodeHg19Male_v0.7.10/encodeHg19Male_bwa-0.7.10.fa
USE_IDR_NBOLEY=true
```

### For cluster use (Sun Grid Engine only)

Add "-s sge" to the command line.

```
$bds -s sge tf_chipseq.bds ...
```

### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
