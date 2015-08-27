ATAC Seq Pipeline
===================================================


### For Kundaje lab clusters

Environment variables will be automatically set in Kundaje lab clusters.

Add '-kundaje_lab true' at the end of the parameters.
```
$bds atac.bds [BOWTIE_INDEX] [READ1] [READ2] [NTHREADS_BWT2] [GENOMESIZE] [CHROMSIZES] [VPLOT_INDEX] [OUTPUT_DIR] -kundaje_lab true
```
or

Add 'KUNDAJE_LAB= true' to a configuration file.
```
$bds atac.bds [CONF_FILE]

$cat [CONF_FILE]
...
KUNDAJE_LAB= true
...
```

### Parameters from command line arguments 1

```
$bds atac.bds [BOWTIE_INDEX] [READ1] [READ2] [NTHREADS_BWT2] [GENOMESIZE] [CHROMSIZES] [VPLOT_INDEX] [OUTPUT_DIR] -mod [MOD_DEF] -shcmd [ADDITIONAL_INIT] -addpath [PATH_FOR_SOFTWARES]

# If you already have -V option (pass all env. vars to qsub) in your ~/.bds/bds.config and defined all env. vars on your current shell, you can skip these additional parameters.
# Otherwise you need to define enviromnet variables with -mod, -shcmd and -addpath.
```

### Detailed help

To get more detailed help, run atac.bds without any parameters.

```
$bds atac.bds
```


### How to set shell environments (What are MOD, SHCMD and ADDPATH?)

It is important to define enviroment variables (like $PATH) to make bioinformatics softwares in the pipeline work properly. MOD, SHCMD and ADDPATH are three convenient ways to define environment variables. Environment variables defined with MOD, SHCMD and ADDPATH are preloaded for all tasks on the pipeline. For example, if you define environment variables for bwa/0.7.3 with MOD. bwa of version 0.7.3 will be used throughout the whole pipeline (including bwa aln, bwa same and bwa sampe).

1) MOD

There are different versions of bioinformatics softwares (eg. samtools, bedtools and bwa) and <a href="http://modules.sourceforge.net/">Enviroment Modules</a> is the best way to manage environemt variables for them. For example, if you want to add environment variables for bwa 0.7.3 by using Environment Modules. You can simply type the following:

```
$ module add bwa/0.7.3;
```

The equivalent setting in the pipeline configuration file should look like:
```
MOD= bwa/0.7.3;
```

You can have multiple lines for MOD since any suffix is allowed. Use ; as a delimiter.
```
MOD_BIO= bwa/0.7.3; bedtools/2.x.x; samtools/1.2
MOD_LANG= r/2.15.1; java/latest
```

2) SHCMD

If you have softwares locally installed on your home, you may need to add to them environment variables like $PATH, $LD_LIBRARY_PATH and so on. <b>IMPORTANT!</b> Note that any pre-defined enviroment variables (like $PATH) should be referred in a curly bracket like ${PATH}. This is because BDS distinguishes environment variables from BDS variables by a curly bracket ${}.
```
SHCMD= export PATH=${PATH}:path_to_your_program
```

You can have multiple lines for SHCMD since any suffix is allowed. Use ; as a delimiter. 
```
SHCMD_R= export PATH=${PATH}:/home/userid/R-2.15.1;
SHCMD_LIB= export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/R-2.15.1/lib
```

SHCMD is not just for adding environemt variables. It can execute any bash shell commands prior to any jobs on the pipeline. For example, to give all jobs peaceful 10 seconds before running.
```
SHCMD_SLEEP_TEN_SECS_FOR_ALL_JOBS= echo "I am sleeping..."; sleep 10
```

3) ADDPATH

If you just want to add something to your $PATH, use ADDPATH instead of SHCMD. It's much simpler. Use : or ; as a delimiter.

```
ADDPATH= ${HOME}/program1/bin:${HOME}/program1/bin:${HOME}/program2/bin:/usr/bin/test
```


### What are -mod, -shcmd and -addpath?

They are command line argument versions of MOD, SHCMD and ADDPATH. For example,

```
$ bds atac.bds -mod 'bwa/0.7.3; samtools/1.2' -shcmd 'export PATH=${PATH}:/home/userid/R-2.15.1' -addpath '${HOME}/program1/bin'
```




### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
