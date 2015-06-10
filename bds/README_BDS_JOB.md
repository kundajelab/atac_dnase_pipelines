BigDataScript (BDS) pipelines
===================================================================

## BigDataScript

Take a look at <a href="./README_BDS.md">./README_BDS.md</a> for installation, configuration and tips for BigDataScript.

### Baseline pipeline (pipeline.bds)

All BDS pipelines on this repository are based on the general baseline pipeline (pipeline.bds). The baseline pipeline does not involve any bioinformatics analysis. It just parses parameters from cmd. line argument or from configuration files. Important parameters like cpu, mem and walltime and initialization settings (env. module init. and exporting env. vars) should be defined here. There are two ways to define such parameters.

1) From command line arguments 

```
$bds pipeline.bds [OPTS_FOR_PIPELINE]

# example cmd. line args

$bds pipeline.bds -prefix TEST -o out -tmp tmp -wt 7200 ...

```

2) From configuration file
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

### Parameters from configuration file

For any parameters not defined in configuration file, default value will be used. Parameters defined in configuration file overrides those defined in cmd. line arguments.

```
	PREFIX 		: Prefix for all outputs.
	OUTPUT_DIR 	: Output directory. (default: out)
	TMP_DIR 	: Temporary directory for intermediate files. (default: tmp).

	WALLTIME 	: Default walltime in seconds for all cluster jobs (default: 36000).
	NTHREADS 	: Default number of threads for all cluster jobs (default: 1).
	MEMORY 		: Default max. memory in MB for all cluster jobs (default: 4000).

	MODULE 		: Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	SHELLCMD	: Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test")
	ADDPATH		: Paths to be added to env. var. PATH separated by ; or :. (a quicker way to add PATH)

```

### Parameters from command line arguments

For any parameters not defined in cmd. line arguments, default value will be used. The following command will show help for all command line arguments.

```
$bds pipeline.bds -h

	-c <string>       : Configuration file path (if not specified, define parameters in command line argument).

	-prefix <string>  : Prefix for all outputs.
	-o <string>       : Output directory. (default: out)
	-tmp <string>     : Temporary directory for intermediate files. (default: tmp).

	-wt <int>         : Default walltime in seconds for all cluster jobs (default: 36000).
	-nth <int>        : Default number of threads for all cluster jobs (default: 1).
	-mem <int>        : Default max. memory in MB for all cluster jobs (default: 4000).

	-mod <string>     : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	-shcmd <string>   : Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test")
	-addpath <string> : Paths to be added to env. var. PATH separated by ; or :. (a quicker way to add PATH)
```

### How to manage BDS jobs?


### Check BDS job status (HTML progress report) on WEB

If you have a NFS mounted web folder like (/srv/www/kundaje/leepc12/job_status/) connected to the cluster, you can automatically make sync BDS job status (HTML,pdf,log,png,qc) to web folders with the following command.

```
./_tools/sync_bds_report.sh [SRC] [DEST]

# or

./_tools/recursive_ln.sh [SRC] [DEST]

```

You can automatically and recursively sync important status files (html,js,pdf,log,qc) on your working folder to web folder by adding the following to crontab -e.

```
# This is an example for nandi cluster. Updates important files every 5 minutes.
# Modify path for sync_bds_report.sh, SRC (top of your working directory) and DEST (web directory)

*/5 * * * * /users/leepc12/code/pipelines/bds/_tools/sync_bds_report.sh /srv/scratch/leepc12/run /srv/www/kundaje/leepc12/bds_monitor/nandi
```

Or, you can automatically and recursively make symlinks for ALL files on your working folder to web folder by adding the following to crontab -e.

```
# This is an example for nandi cluster. Updates symlinks every 3 minutes.
*/3 * * * * /users/leepc12/code/pipelines/bds/_tools/recursive_ln.sh /srv/scratch/leepc12/run /srv/www/kundaje/leepc12/bds_monitor/nandi
```

sync_bds_report.sh recursively loops through all important status files, make the same directory structure on [DEST] and then copy files to [DEST].
recursive_ln.sh recursively loops through all files,make the same directory structure on [DEST] and then make symlinks on [DEST].

### How to specify softwares for the pipeline?

With a configuration file, you can specify softwares with MODULE*, SHELLCMD* and ADDPATH*.
Any suffix is allowed like MODULE_BIO, MODULE_LANG, SHELLCMD_TOOL, SHELLCMD_TEST and ADDPATH_MATHLIB.

For example, to use bwa 0.7.7 and bedtools 2.x.x and samtools 1.2

```
MODULE_BWA= bwa/0.7.7 
MODULE_ANY_SUFFIX_SHOULD_BE_OKAY= bedtools/2.x.x; samtools/1.2
```

Additional environment variables can be defined with SHELLCMD*

```
# YOU NEED TO COVER ANY ENVIRONMENT VARS WITH CURLY BRACKETS ${} !!
SHELLCMD_ETC= export PATH="${PATH}:/usr/bin/example"
SHELLCMD_TEST= TEST_PATH="/usr/lib/example"; export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${TEST_PATH}"
SHELLCMD_YOU_CAN_DEFINE_ANY_CUSTOM_VAR= CUST_VAR=200
```

With cmd. line arguments, use -mod [], -shcmd [] and -addpath []. Use semicolon as a delimiter.

```
bds $BDS_SCRIPT -mod 'bwa/0.7.7; bedtools/2.x.x; samtools/1.2' -shcmd 'export PATH="${PATH}:/usr/bin/example"; ...'
```

Running BDS script with the above configuration file or cmd. line arguments will actually execute the following:

```
module add bwa/0.7.7
module add bedtools/2.x.x
module add samtools/1.2

export PATH="${PATH}:/usr/bin/example"
TEST_PATH="/usr/lib/example"
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${TEST_PATH}"
CUST_VAR=200
```


### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
