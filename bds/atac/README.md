ATAC Seq Pipeline
===================================================


### For Kundaje lab memebers

Environment variables will be automatically set in Kundaje lab clusters.

Add '-kundaje_lab true' at the end of the parameters.
```
$bds atac.bds [BOWTIE_IDX] [READ1] [READ2] [NUMTHREADS] [GENOMESIZE] [CHROMSIZE] [V_INDEX] [OUTPUTDIR] -kundaje_lab true
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
$bds atac.bds [BOWTIE_IDX] [READ1] [READ2] [NUMTHREADS] [GENOMESIZE] [CHROMSIZE] [V_INDEX] [OUTPUTDIR] -mod "${MOD_DEF}" -shcmd "{ADDITIONAL_INIT}" -addpath "{PATH_FOR_SOFTWARES}"

# If you already have -V option (pass all env. vars to qsub) in your ~/.bds/bds.config and defined all env. vars on your current shell, you can skip these additional parameters.
# Otherwise you need to define enviromnet variables with -mod, -shcmd and -addpath.

# For example (IMPORTANT! use -mod "${MOD_DEF}" instead of -mod ${MOD_DEF}. Without "", each blank will play as a delimiter so will get errors.

MOD_DEF= 'bowtie/2.2.4; samtools/1.2; bedtools/2.21.0; picard-tools/1.129; ucsc_tools/3.0.9; MACS2/2.1.0; java/latest; preseq/1.0.2'
PATH_FOR_SOFTWARES= '/biotools/bin:/whatever/software/you/need'
ADDITIONAL_INIT= 'export _JAVA_OPTIONS="-Xms256M -Xmx512M -XX:ParallelGCThreads=1"; export MAX_JAVA_MEM="4G"; export MALLOC_ARENA_MAX=4'

# There is another nice example on ./conf/cmds_all_jobs_atac_eugene_CORRECTED.sh, which includes more than 20 samples.
```

### Parameters from command line arguments 2

```
$bds atac.bds [OPTS_FOR_ATAC]

#[OPTS_FOR_ATAC] are like the following:

	-KUNDAJE_LAB <bool>   : Set it true if you run the pipeline on Kundaje lab servers (automatically set environments, default: false)

	-c <string>           : Configuration file path (if not specified, define parameters in command line argument).

	-prefix <string>      : Prefix for all outputs.
	-o <string>           : Output directory. (default: out)
	-tmp <string>         : Temporary directory for intermediate files. (default: tmp).

	-gen <string>         : Reference genome name for epigenome browser track generation (eg. hg19, hg18, mm9 or mm10)

	-wt <int>             : Default walltime in seconds for all cluster jobs (default: 36000).
	-nth <int>            : Default number of threads for all cluster jobs (default: 1).
	-mem <int>            : Default max. memory in MB for all cluster jobs (default: 4000).

	-mod <string>         : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	-shcmd <string>       : Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test").
	-addpath <string>     : Paths to be added to env. var. PATH separated by ; or :. (a quicker way to add PATH)

	-BOWTIE_IDX <string>  : Path for bowtie index
	-READ1 <string>       : Read1 fastq
	-READ2 <string>       : Read2 fastq
	-NUMTHREADS <int>     : Number of threads for bowtie2
	-GENOMESIZE <string>  : 'hs' by default
	-CHROMSIZE <string>   : Path for chrom.sizes file for your .fa
	-V_INDEX <string>     : Index for v-plot
```

### Parameters from configuration file

```
$bds atac.bds [CONF_FILE]

$cat [CONF_FILE]

	KUNDAJE_LAB : Set it true if you run the pipeline on Kundaje lab servers (automatically set environments, default: false)

	PREFIX      : Prefix for all outputs.
	OUTPUT_DIR  : Output directory. (default: out)
	TMP_DIR     : Temporary directory for intermediate files. (default: tmp).

	REF_GENOME  : Reference genome name for epigenome browser track generation (eg. hg19, hg18, mm9 or mm10)

	WALLTIME    : Default walltime in seconds for all cluster jobs (default: 36000).
	NTHREADS    : Default number of threads for all cluster jobs (default: 1).
	MEMORY      : Default max. memory in MB for all cluster jobs (default: 4000).

	MODULE      : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	SHELLCMD    : Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test")
	ADDPATH     : Paths to be added to env. var. PATH separated by ; or :. (a quicker way to add PATH)

	BOWTIE_IDX  : Path for bowtie index
	READ1       : Read1 fastq
	READ2       : Read2 fastq
	NUMTHREADS  : Number of threads for bowtie2
	GENOMESIZE  : 'hs' by default
	CHROMSIZE   : Path for chrom.sizes file for your .fa
	V_INDEX     : Index for v-plot
```


### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
