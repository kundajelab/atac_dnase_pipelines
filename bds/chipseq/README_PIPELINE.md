BigDataScript (BDS) Pipelines
===============================================




### Installation instruction

Get BigDataScript (v0.9999 is stable and doesn't require high java version).
```
$ git clone https://github.com/pcingola/BigDataScript
$ cd BigDataScript
$ git checkout tags/v0.9999
$ cp distro/bds_Linux.tgz $HOME
$ cd $HOME
$ tar zxvf bds_Linux.tgz
$ rm bds_Linux.tgz
```

Find bds.config and move it to $HOME/.bds/.
```
$ mkdir -p $HOME/.bds
$ cp bds.config $HOME/.bds/
```

(Recommended if you don't know what these are) Add the following lines to your $HOME/.bashrc or $HOME/.bash_profile:
```
export _JAVA_OPTIONS="-Xms256M -Xmx512M -XX:ParallelGCThreads=1"
export MAX_JAVA_MEM="16G"
export MALLOC_ARENA_MAX=4
export PATH=$PATH:$HOME/.bds
```



### Installation on SCG3 (scg3, carmack and crick)

All dependencies except BDS have already been installed on SCG3. Follow the instruction in the previous section.




### Installation on Kundaje lab clusters (mitra, nandi, vayu, kali, durga, wotan and amold)

All dependencies have already been installed on lab servers. Find bds.config in the pipeline repo and move it to $HOME/.bds/.
```
$ mkdir -p $HOME/.bds
$ cp bds.config $HOME/.bds/
```



### How to run pipelines?

Create and move to your working directory. Run the following command:
```
$ bds [PIPELINE_BDS] [PARAMETERS]
```

Do not run multiple BDS pipelines on the same working directory. BDS creates an HTML report and temporary files on the working directory. Things will be messed up.




### Running pipelines with a cluster engine

You can run BDS pipeline with a specified cluster engine. Choose your cluster system (local: UNIX threads, sge: Sun Grid Engine, ...).
```
$ bds -s [SYSTEM] [PIPELINE.BDS] ...
```

Modify `$HOME./.bds/bds.config` to change your default system. The following example is to use Sun Grid Engine (sge) as your default system. Then you no longer need to add `-s sge` to the command line.
```
#system = local
system = sge
```

You need additional modification on bds.config to correctly configure your cluster engine. Read more on <a href="http://pcingola.github.io/BigDataScript/bigDataScript_manual.html" target="_blank">http://pcingola.github.io/BigDataScript/bigDataScript_manual.html</a>. For Kundaje lab clusters and SCG3, it's already set up for Sun Grid Engine.




### How to display all parameters and help

Run pipelines without parameters.
```
$ bds [PIPELINE_BDS]
```





### Resource settings (walltime and max. memory)

Most clusters have resource limitation so that jobs submitted without it will be declined. By default, walltime is 11 hours and max memory is 8GB. To change them, add the following parameters to the command line.
```
-wt [WALLTIME; examples: 13:20:20, 10h, 7200] -mem [MAX_MEMORY; examples: 5G, 2000K]
```

You can also specify walltime and max. memory for a specific job. To see which job has specific resource settings, run the pipeline without parameters `$ bds [PIPELINE_BDS]` then it will display all parameters including resource settings and help. The following line is an example parameter to increase walltime and max. memory for MACS2 peak calling.
```
-wt_macs2 40h -mem 20G
```

If your system (either local or cluster engine) doesn't limit walltime and max. memory for jobs, add the following to the command line. Pipeline jobs will run without resource restriction.
```
-skip_wt_mem
```



### Resource settings On SCG3

You always need to submit pipeline jobs to Sun Grid Engine. Carefully define resources settings on SCG3. If walltime is over 6 hours on SCG3, jobs will be submitted to a longer queue, which makes you wait long to get jobs executed. By default, walltime is 11 hours. If your job is small enough to be done in 6 hours, then make the walltime shorter than 6h.
```
-wt 5:50:00
```



### Resource settings on Kundaje lab clusters

There is no limit for walltime and max. memory on Kundaje lab clusters.





### How to stop pipelines?

A Pipeline will go through several stages. To stop it, just press `Ctrl+C`. It will stop the whole pipeline and remove all intermediate files in its output directory. It doesn't remove outputs from stages which are succesfully finished. It doesn't affect jobs on other BDS pipelines. 

You can resume the pipeline by using the same command line you used to start it. For each stage, BDS automatically check if output files exist or they are newer than input files, and then determine whether stages need to be re-run or not. 





### How to deal with BDS pipeline errors?


1) Take a look at HTML report (which contains all STDERR/STDOUT for all jobs in the pipeline). It tells you everything about all pipeline jobs. Find which stage is errorneous. Carefully look at system messages (STDERR and STDOUT) for it.

2) Correct errors.
   2-1) Lack of memory: increase memory for all jobs (e.g. add -mem 20G) or a specific problematic job (e.g. add -mem_macs2 20G).
   2-2) Timeout: increase walltime for all jobs (e.g. add -wt 24h) or a specific long job (e.g. add -wt_macs2 200h).
                 (Warning! Most clusters have limit for walltime. Make it as shortest as you can to get your queued jobs executed quickly.)
   2-3) Wrong input: check all input files are available.
   2-4) Software error: use recommended software versions.

3) Resume pipeline with the same command line that you used to start it with. Previous successful stages will be automatically skipped.





### How to define parameters?

There are two ways to define parameters for pipelines. Default values are already given for most of parameters. Take a look at example commands and configuration files (./examples). Both methods share the same key names.

1) From command line arguments 
```
$ bds [PIPELINE_BDS] [PARAMETERS]
```
Example for chipseq pipeline:
```
$ bds chipseq.bds \
-fastq1 /DATA/ENCFF000YLW.fastq.gz \
-fastq2 /DATA/ENCFF000YLY.fastq.gz \
-ctl_fastq1 /DATA/Ctl/ENCFF000YRB.fastq.gz \
-bwa_idx /INDEX/encodeHg19Male_bwa-0.7.3.fa
...
```

2) From a configuration file
```
$ bds [PIPELINE_BDS] [CONF_FILE]
```
Example configuriation file for chipseq pipeline:
```
$ cat [CONF_FILE]

fastq1= /DATA/ENCFF000YLW.fastq.gz
fastq2= /DATA/ENCFF000YLY.fastq.gz
ctl_fastq1= /DATA/Ctl/ENCFF000YRB.fastq.gz
bwa_idx= /INDEX/encodeHg19Male_bwa-0.7.3.fa
...
```

You can override any parameters defined in a configuration file by adding them to the command line.
```
$ bds [PIPELINE_BDS] [CONF_FILE] [PARAMS_TO_BE_OVERRIDEN]
```



### Using species file

There are many species specific parameters like indices (bwa, bowtie, ...), chromosome sizes and sequence files (chr*.fa). If you have multiple pipelines, it's inconvenient to individually define all parameters in a command line argument for each pipeline run. However, if you have a species file with all species specific parameters defined, then you define less parameters in the command line and share the species file with all other pipelines.

Add the following to the command line to specify species and species file.
```
-species [SPECIES; hg19, mm9, ...] -species_file [PATH_FOR_SPECIES_FILE]
```

You can override any parameters defined in the species file by adding them to command line argument or configuration file. For example, if you want to override parameters for BWA index and umap:
```
-species hg19 -species_file my_species.conf -bwa_idx [YOUR_OWN_BWA_IDX] -chrsz [YOUR_OWN_CHR_SIZES_FILE]
```

Example species file looks like the following. You can define your own species.
```
[hg19]
chrsz   = /mnt/data/annotations/by_release/hg19.GRCh37/hg19.chrom.sizes // chrome sizes
seq     = /mnt/data/ENCODE/sequence/encodeHg19Male 			// genome reference sequence
gensz   = hs // genome size: hs for humna, mm for mouse
umap    = /mnt/data/ENCODE/umap/encodeHg19Male/globalmap_k1tok1000 	// uniq. mappability tracks
umap_hic= /mnt/data/ENCODE/umap/encodeHg19Male/globalmap_k20tok54 	// uniq. mappability tracks
bwa_idx = /mnt/data/annotations/indexes/bwa_indexes/encodeHg19Male/v0.7.10/encodeHg19Male_bwa-0.7.10.fa
bwt_idx = /mnt/data/annotations/indexes/bowtie1_indexes/encodeHg19Male/encodeHg19Male
bwt2_idx = /mnt/data/annotations/indexes/bowtie2_indexes/bowtie2/ENCODEHg19_male
vplot_idx = /mnt/data/annotations/indexes/vplot_indexes/hg19/parsed_hg19_RefSeq.merged.ANS.bed

# your own definition for species
[hg19_custom]
chrsz   = ...
seq     = ...
...

[mm9]
...

[mm10]
...
```

Description for parameters in a species file.
```
chrsz               : Chrome sizes file path (use fetchChromSizes from UCSC tools).
seq                 : Reference genome sequence directory path (where chr*.fa exist).
gensz               : Genome size; hs for human, mm for mouse (default: hs).
umap                : Unique mappability tracks directory path (https://sites.google.com/site/anshulkundaje/projects/mappability).
umap_hic            : Unique mappability tracks directory path (for HiC, DO NOT USE all-mappable umap track starting from 1bp)
bwa_idx             : BWA index (full path prefix of [].bwt file) .
bwt_idx             : Bowtie index (full path prefix of [].1.ebwt file).
bwt2_idx            : Bowtie2 index (full path prefix of [].1.bt2 file).
vplot_idx           : V plot index (full path for bed file).
```




### Using species file on SCG3 and Kundaje lab clusters

Add the following to the command line and that's it.
```
-species [SPECIES; hg19, mm9]
```

hg19 and mm9 are available for SCG3 and Kundaje lab clusters. If you are interested in other species, add species to `species/species*.conf` and share with lab members or create your own species file.





### Useful HTML report for debugging

BDS HTML report is located at the working folder with name [PIPELINE_NAME]_[TIMESTAMP]_report.html. This report is automatically generated by BDS and is very useful for debugging since it shows summary, timeline, Stdout and Stderr for all jobs in the pipeline.





### Dry run

Dry run (this actually does nothing) to check next stages and input/output file names for it.
```
$ bds -dryRun [PIPELINE_BDS] ...
```




### How to set shell environments

Ignore this section if you are working on SCG3 or Kundaje lab clusters.

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
shcmr_R= export PATH=${PATH}:/home/userid/R-2.15.1;
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
$ bds [PIPELINE_BDS] -mod 'bwa/0.7.3; samtools/1.2' -shcmd 'export PATH=${PATH}:/home/userid/R-2.15.1' -addpath '${HOME}/program1/bin'
```



### Troubleshooting

If see the following error when you submit jobs to Sun Grid Enginee,
```
/bin/bash: module: line 1: syntax error: unexpected end of file
```

Check your $HOME/.bashrc if it has any errorneous lines.

Remove the following line in you module initialization scripts ($modULESHOME/init/bash or /etc/profile.d/modules.sh).
```
export -f module
```




### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
