ATAC Seq Pipeline
===================================================

### Installation instruction

Get BigDataScript (v0.9999).
```
$ git clone https://github.com/pcingola/BigDataScript
$ cd BigDataScript
$ git checkout tags/v0.9999
$ cp distro/bds_Linux.tgz $HOME
$ cd $HOME
$ tar zxvf bds_Linux.tgz
```

Get the latest version of ATAC pipeline.
```
$ git clone https://github.com/kundajelab/pipelines/
$ cd bds/atac
$ mkdir -p $HOME/.bds
$ cp bds.config $HOME/.bds/
```

Add the following lines to your $HOME/.bashrc or $HOME/.bash_profile:
```
export _JAVA_OPTIONS="-Xms256M -Xmx512M -XX:ParallelGCThreads=1"
export MAX_JAVA_MEM="8G"
export MALLOC_ARENA_MAX=4
export PATH=$PATH:$HOME/.bds
```

### Installation instruction (for Kundaje lab members)

For Kundaje lab members, BDS and all dependencies have already been installed on lab servers. Do not run install_dependencies.sh on Kundaje lab servers.
```
$ mkdir -p ~/.bds
$ cp bds.config ~/.bds
```

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


### Using Species file

For the pipeline, there are many species specific parameters like indices (bwa, bowtie, ...), chrome sizes, sequence file and genome size. If you have multiple pipelines, it's a hard job to individually define all parameters for each pipeline. However, if you have a species file with all species specific parameters defined, then you define less parameters and share the species file with all other pipelines.

If species file is not defined, pipeline looks for paths in the following order:
```
1) Configruation file for the pipeline
2) Species file defined by '-species_file' option
3) [BDS_SCRIPT_PATH]/species.conf
4) [WORK_DIR]/species.conf
```
If -kundaje_lab is define.
```
5) [BDS_SCRIPT_PATH]/species/species_kundaje_lab.conf
```

You can override any parameters defined in the species file by adding them to command line argument or configuration file.
```
$ bds chipseq.bds ... -species [SPECIES] -species_file [SPECIES_FILE] ... [ANY_PARAMETETRS_TO_BE_OVERRIDEN]
```

For example, if you want to override parameters for BWA index and umap:
```
$ bds chipseq.bds ... -species [SPECIES] -species_file [SPECIES_FILE] ... -bwa_idx [YOUR_OWN_BWA_IDX] -umap [YOUR_OWN_UMAP]
```

Example species file:
```
[hg19]
chrsz   = /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes // chrome sizes
seq     = /mnt/data/ENCODE/sequence/encodeHg19Male // genome reference sequence
gensz   = hs // genome size: hs for humna, mm for mouse
umap    = /mnt/data/ENCODE/umap/encodeHg19Male/globalmap_k20tok54 // uniq. mappability tracks
bwa_idx = /mnt/data/annotations/indexes/bwa_indexes/encodeHg19Male/v0.7.10/encodeHg19Male_bwa-0.7.10.fa
bwt_idx = /mnt/data/annotations/indexes/bowtie1_indexes/encodeHg19Male/encodeHg19Male
bwt2_idx = /mnt/data/annotations/indexes/bowtie2_indexes/bowtie2/ENCODEHg19_male
vplot_idx = /mnt/data/annotations/indexes/vplot_indexes/hg19/parsed_hg19_RefSeq.merged.ANS.bed

[hg38]
...

[mm9]
...

[mm10]
...
```


### Detailed help

To get more detailed help, run atac.bds without any parameters.

```
$ bds atac.bds
```

### For cluster use (Sun Grid Engine only)

Add "-s sge" to the command line.

```
$ bds -s sge atac.bds [...]
```

### How to set shell environments (What are mod, shcmd and addpath?)

It is important to define enviroment variables (like $PATH) to make bioinformatics softwares in the pipeline work properly. mod, shcmd and addpath are three convenient ways to define environment variables. Environment variables defined with mod, shcmd and addpath are preloaded for all tasks on the pipeline. For example, if you define environment variables for bwa/0.7.3 with mod. bwa of version 0.7.3 will be used throughout the whole pipeline (including bwa aln, bwa same and bwa sampe).

1) mod

There are different versions of bioinformatics softwares (eg. samtools, bedtools and bwa) and <a href="http://modules.sourceforge.net/">Enviroment Modules</a> is the best way to manage environemt variables for them. For example, if you want to add environment variables for bwa 0.7.3 by using Environment Modules. You can simply type the following:

```
$ module add bwa/0.7.3;
```

The equivalent setting in the pipeline configuration file should look like:
```
mod= bwa/0.7.3;
```

You can have multiple lines for mod since any suffix is allowed. Use ; as a delimiter.
```
mod_bio= bwa/0.7.3; bedtools/2.x.x; samtools/1.2
mod_lang= r/2.15.1; java/latest
```

2) shcmd

If you have softwares locally installed on your home, you may need to add to them environment variables like $PATH, $LD_LIBRARY_PATH and so on. <b>IMPORTANT!</b> Note that any pre-defined enviroment variables (like $PATH) should be referred in a curly bracket like ${PATH}. This is because BDS distinguishes environment variables from BDS variables by a curly bracket ${}.
```
shcmd= export PATH=${PATH}:path_to_your_program
```

You can have multiple lines for shcmd since any suffix is allowed. Use ; as a delimiter. 
```
shcmd_r= export PATH=${PATH}:/home/userid/R-2.15.1;
shcmd_lib= export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/R-2.15.1/lib
```

shcmd is not just for adding environemt variables. It can execute any bash shell commands prior to any jobs on the pipeline. For example, to give all jobs peaceful 10 seconds before running.
```
shcmd_SLEEP_TEN_SECS_FOR_ALL_JOBS= echo "I am sleeping..."; sleep 10
```

3) addpath

If you just want to add something to your $PATH, use addpath instead of shcmd. It's much simpler. Use : or ; as a delimiter.

```
addpath= ${HOME}/program1/bin:${HOME}/program1/bin:${HOME}/program2/bin:/usr/bin/test
```


### What are -mod, -shcmd and -addpath?

They are command line argument versions of mod, shcmd and addpath. For example,

```
$ bds atac.bds -mod 'bwa/0.7.3; samtools/1.2' -shcmd 'export PATH=${PATH}:/home/userid/R-2.15.1' -addpath '${HOME}/program1/bin'
```




### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
