HiC Pipeline
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

Get the latest version of HiC pipeline.
```
$ git clone https://github.com/kundajelab/pipelines/
$ cd bds/hic
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

For Kundaje lab members, add the following to the command line instead of defining all species-specific parameters like BWA index and umap path.
```
-kundaje_lab -species [SPECIES: hg19, mm9, ...]
```

### How does HiC pipeline works?

There are two stages for HiC pipeline.

1) Mapping and sorting (hic_map.bds)
: Inputs are fastqs. Map, align and sort them to generate cleaned pairs.

2) HiC (hic.bds)
: Using cleaned pairs from the stage 1), perform HiC analysis.

Two stages share the same indices for librarie (lib) and replicates (rep).


### Usage (mapping and sorting)

The following examples are for libraries (libs: D3, D6 and SC) and replicates (reps: R2 and R8).

1) Define parameters in command line argument. 
```
$ bds hic_map.bds \
-libs D3,D6,SC \
-reps R2,R8 \
-fastq_L[Lib_ID]_R[Rep_ID]_P[Pair_ID] [FASTQ_FOR_LIB1_REP1_PAIR1] \
...
-bwa_idx [BWA_INDEX]
-nth_bwa [# THREADS FOR BWA]
```

2) Define parameters in configuration file.
```
$ bds hic_map.bds [CONF_FILE]

$ cat [CONF_FILE]
libs = D3,D6,SC
reps = R2,R8
fastq_L[Lib_ID]_R[Rep_ID]_P[Pair_ID] = [FASTQ_FOR_LIB1_REP1_PAIR1]
...
bwa_idx = [BWA_INDEX]
nth_bwa = [# THREADS FOR BWA]
```


### Usage (HiC)

1) Using Root directory of mapped data
```
$ bds hic.bds \
-libs D3,D6,SC \
-reps R2,R8 \
-root_mapped [ROOT_OF_MAPPED_DATA: OUT_DIR OF MAPPING STAGE] \
...
-merge [MERGE_METHOD; 0:no_merge, 1:merge_replicates_only, 2:merge_libraries_and_replicates] \
-res [COMMA-SEPERATED RESOLUTIONS; Example: 100,200,1000] \
-blacklist [BLACKLIST_FILE] \
-RE_file [RE_FILE] \
-umap [MAPPABILITY_DATA]
```

2) Using paths for individual cleaned pairs
```
$ bds hic.bds \
-libs D3,D6,SC \
-reps R2,R8 \
-cln_pair_L1_R1 [CLEANED_PAIR_FOR_LIB1_REP1] \
...
-merge [MERGE_METHOD; 0:no_merge, 1:merge_replicates_only, 2:merge_libraries_and_replicates] \
-res [COMMA-SEPERATED RESOLUTIONS; Example: 100,200,1000] \
-blacklist [BLACKLIST_FILE] \
-RE_file [RE_FILE] \
-umap [MAPPABILITY_DATA]
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
If -kundaje_lab is defined.
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

To get more detailed help, run hic.bds without any parameters.

```
$ bds hic_map.bds

$ bds hic.bds
```

### For cluster use (Sun Grid Engine only)

Add "-s sge" to the command line.

```
$ bds -s sge hic_map.bds [...]

$ bds -s sge hic.bds [...]
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
$ bds hic.bds -mod 'bwa/0.7.3; samtools/1.2' -shcmd 'export PATH=${PATH}:/home/userid/R-2.15.1' -addpath '${HOME}/program1/bin'
```




### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
