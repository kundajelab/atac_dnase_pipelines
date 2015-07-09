BigDataScript (BDS) pipelines
===================================================================

## BigDataScript

Take a look at <a href="README_BDS.md">README_BDS.md</a> for installation, configuration and tips (logging/debugging and reporting) for BigDataScript. IMPORTANT! BDS pipelines will not work if you don't have properly installed BDS and its config (bds.config) on your home ($HOME/.bds/bds.config).

There are three bioinformatics pipelines available on this repository.

1) <a href="atac/atac.bds">ATAC-seq</a>
2) <a href="chipseq/tf_chipseq.bds">TF ChIP-seq and Histone ChIP-seq</a>
2) <a href="chipseq/hist_chipseq.bds">TF ChIP-seq and Histone ChIP-seq</a>


### Baseline pipeline (pipeline.bds)

All BDS pipelines on this repository are based on the general baseline pipeline (pipeline.bds). The baseline pipeline does not involve any bioinformatics analysis. It just parses parameters from command line argument or from a configuration file. Important parameters like cpu, mem and walltime and initialization settings (env. module init. and exporting env. vars) should be defined here. There are two ways to define such parameters.

1) From command line arguments 
```
$bds pipeline.bds [OPTS_FOR_PIPELINE]

# example cmd. line args

$bds pipeline.bds -prefix TEST -o out -tmp tmp -wt 7200 ...

```

2) From a configuration file
```
$bds pipeline.bds -c [CONF_FILE]

# or

$bds pipeline.bds [CONF_FILE]

# example configuration file
$cat [CONF_FILE]

PREFIX= TEST
OUTPUT_DIR= out
TMP_DIR= tmp
WALLTIME= 7200 
...
```


### Parameters from a configuration file

For any parameters not defined in a configuration file, default value will be used. Parameters defined in a configuration file overrides those defined in command line arguments.

```
	PREFIX 		: Prefix for all outputs.
	OUTPUT_DIR 	: Output directory. (default: out)
	TMP_DIR 	: Temporary directory for intermediate files. (default: tmp).

	REF_GENOME	: Reference genome name for epigenome browser track generation (eg. hg19, hg18, mm9 or mm10)

	WALLTIME 	: Default walltime in seconds for all cluster jobs (default: 36000).
	NTHREADS 	: Default number of threads for all cluster jobs (default: 1).
	MEMORY 		: Default max. memory in MB for all cluster jobs (default: 4000).

	MODULE 		: Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	SHELLCMD	: Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test")
	ADDPATH		: Paths to be added to env. var. PATH separated by ; or :. (a quicker way to add PATH)

```


### Parameters from command line arguments

For any parameters not defined in command line arguments, default value will be used. The following command will show help for all command line arguments.

```
$bds pipeline.bds -h

        -c <string>        : Configuration file path (if not specified, define parameters in command line argument).

        -prefix <string>   : Prefix for all outputs.
        -o <string>        : Output directory. (default: out)
        -tmp <string>      : Temporary directory for intermediate files. (default: tmp).

        -gen <string>      : Reference genome name for epigenome browser track generation (eg. hg19, hg18, mm9 or mm10)

        -wt <int>          : Default walltime in seconds for all cluster jobs (default: 22000).
        -nth <int>         : Default number of threads for all cluster jobs (default: 1).
        -mem <int>         : Default max. memory in MB for all cluster jobs (default: 8000).

        -mod <string>      : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
        -shcmd <string>    : Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test").
        -addpath <string>  : Path to be appended to env. var. PATH. Multiple paths should be separated by ; (example: "/bin/test:/bin/test2")
```


## What are MODULE, SHELLCMD and ADDPATH?

It is important to define enviroment variables (like $PATH) to make bioinformatics softwares in the pipeline work properly. MODULE, SHELLCMD and ADDPATH are three convenient ways to define environment variables. Environment variables defined with MODULE, SHELLCMD and ADDPATH are preloaded for all tasks on the pipeline. For example, if you define environment variables for bwa/0.7.10 with MODULE. bwa of version 0.7.10 will be used throughout the whole pipeline (including bwa aln, bwa same and bwa sampe).

1) MODULE

There are different versions of bioinformatics softwares (eg. samtools, bedtools and bwa) and <a href="http://modules.sourceforge.net/">Enviroment Modules</a> is the best way to manage environemt variables for them. For example, if you want to add environment variables for bwa 0.7.10 by using Environment Modules. You can simply type the following:

```
$module add bwa/0.7.10;
```

The equivalent setting in the pipeline configuration file should look like:
```
MODULE= bwa/0.7.10;
```

You can have multiple lines for MODULE since any suffix is allowed. Use ; as a delimiter.
```
MODULE_BIO= bwa/0.7.10; bedtools/2.x.x; samtools/1.2
MODULE_LANG= r/2.15.1; java/latest
```

2) SHELLCMD

If you have softwares locally installed on your home, you may need to add to them environment variables like $PATH, $LD_LIBRARY_PATH and so on. <b>IMPORTANT!</b> Note that any pre-defined enviroment variables (like $PATH) should be referred in a curly bracket like ${PATH}. This is because BDS distinguishes environment variables from BDS variables by a curly bracket ${}.
```
SHELLCMD= export PATH=${PATH}:path_to_your_program
```

You can have multiple lines for SHELLCMD since any suffix is allowed. Use ; as a delimiter. 
```
SHELLCMD_R= export PATH=${PATH}:/home/userid/R-2.15.1;
SHELLCMD_LIB= export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/R-2.15.1/lib
```

SHELLCMD is not just for adding environemt variables. It can execute any bash shell commands prior to any jobs on the pipeline. For example, to give all jobs peaceful 10 seconds before running.
```
SHELLCMD_SLEEP_TEN_SECS_FOR_ALL_JOBS= echo "I am sleeping..."; sleep 10
```

3) ADDPATH

If you just want to add something to your $PATH, use ADDPATH instead of SHELLCMD. It's much simpler. Use : or ; as a delimiter.

```
ADDPATH= ${HOME}/program1/bin:${HOME}/program1/bin:${HOME}/program2/bin:/usr/bin/test
```

### What are -mod, -shcmd and -addpath?

They are command line argument versions of MODULE, SHELLCMD and ADDPATH. For example,

```
$bds pipeline.bds -mod 'bwa/0.7.10; samtools/1.2' -shcmd 'export PATH=${PATH}:/home/userid/R-2.15.1' -addpath '${HOME}/program1/bin'
```

### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
