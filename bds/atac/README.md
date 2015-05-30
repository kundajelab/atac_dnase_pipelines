ATAC Seq Pipeline
===================================================

Please take a look at ../README.md first.


### Parameters from command line arguments 1

```
bds atac.bds [BOWTIE_IDX] [READ1] [READ2] [NUMTHREADS] [GENOMESIZE] [CHROMSIZE] [OUTPUTDIR]
```

### Parameters from command line arguments 2

```
	-c <string>           : Configuration file path (if not specified, define parameters in command line argument).

	-prefix <string>      : Prefix for all outputs.
	-o <string>           : Output directory. (default: out)
	-tmp <string>         : Temporary directory for intermediate files. (default: tmp).

	-wt <int>             : Default walltime in seconds for all cluster jobs (default: 36000).
	-nth <int>            : Default number of threads for all cluster jobs (default: 1).
	-mem <int>            : Default max. memory in MB for all cluster jobs (default: 4000).

	-mod <string>         : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	-shcmd <string>       : Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test").

	-BOWTIE_IDX <string>  : Path for bowtie index
	-READ1 <string>       : Read1 fastq
	-READ2 <string>       : Read2 fastq
	-NUMTHREADS <int>     : Number of threads for bowtie2
	-GENOMESIZE <string>  : 'hs' by default
	-CHROMSIZE <string>   : Path for chrom.sizes file for your .fa
```

### Parameters from configuration file

```
	PREFIX      : Prefix for all outputs.
	OUTPUT_DIR  : Output directory. (default: out)
	TMP_DIR     : Temporary directory for intermediate files. (default: tmp).

	WALLTIME    : Default walltime in seconds for all cluster jobs (default: 36000).
	NTHREADS    : Default number of threads for all cluster jobs (default: 1).
	MEMORY      : Default max. memory in MB for all cluster jobs (default: 4000).

	MODULE      : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	SHELLCMD    : Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test")

	BOWTIE_IDX  : Path for bowtie index
	READ1       : Read1 fastq
	READ2       : Read2 fastq
	NUMTHREADS  : Number of threads for bowtie2
	GENOMESIZE  : 'hs' by default
	CHROMSIZE   : Path for chrom.sizes file for your .fa
```


### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
