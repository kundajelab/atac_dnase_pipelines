BigDataScript (BDS) pipelines
===================================================================

### Installation instruction for BDS

Check if BDS is already installed on a cluster.

```
which bds
```
If not, install it locally on your home.

```
cd $HOME
wget -O bds_Linux.tgz https://github.com/pcingola/BigDataScript/releases/download/v0.999k/bds_Linux.tgz
tar zxvf bds_Linux.tgz
rm bds_Linux.tgz
```

Download bds.config.SERVERNAME on the git repo. Move it to your $HOME/.bds/bds.config. In bds.config, you can specify your information like email address to get job status notification from clusters.

```
clusterRunAdditionalArgs = -A accountID -M user@gmail.com
```

If you want to pass all of environment variables on your current login shell to cluster job nodes, add -V to clusterRunAdditionalArgs

```
clusterRunAdditionalArgs = -V -A accountID -M user@gmail.com
```

<b> Skip steps below if BDS is already installed on /bin or /usr/bin </b>

Add the following to $HOME/.bds/bds.config. 

```
# This is example, Do not use $HOME in bds.config. Write absolute path for your $HOME/.bds.
clusterRunAdditionalArgs = -v PATH=/users/leepc12/.bds
```

Also, modify your $HOME/.bashrc to add the following:

```
export PATH=$PATH:$HOME/.bds/
```

Official homepage for BDS is at <a href="http://pcingola.github.io/BigDataScript/download.html">http://pcingola.github.io/BigDataScript/download.html</a> and the git repo for BDS is at <a href="https://github.com/pcingola/BigDataScript.git">https://github.com/pcingola/BigDataScript.git</a>.


### Running BDS Script

Running locally,

```
bds $BDS_SCRIPT 
```

Running on a cluster with grid engine.

```
bds -s sge $BDS_SCRIPT 
```

Debugging

```
bds -d $BDS_SCRIPT
```

<b> IMPORTANT </b>

You can give cmd. line argument to $BDS_SCRIPT (not BDS itself)

```
bds $BDS_SCRIPT arg0 arg1 arg2 arg3 ...
```

For example, if you want to run tf_chipseq BDS script on a grid engine.

```
bds -s sge $BDS_SCRIPT tf_chipse_SE_XXXX.conf
```

### Parameters from configuration file

For any parameters not defined in configuration file, default value will be used. Parameters defined in configuration file overrides those defined in cmd. line arguments.

```
	CONF_FILE 	: Configuration file path (if not specified, define parameters in command line argument).

	PREFIX 		: Prefix for all outputs.
	OUTPUT_DIR 	: Output directory. (default: out)
	TMP_DIR 	: Temporary directory for intermediate files. (default: tmp).

	WALLTIME 	: Default walltime in seconds for all cluster jobs (default: 36000).
	NTHREADS 	: Default number of threads for all cluster jobs (default: 1).
	MEMORY 		: Default max. memory in MB for all cluster jobs (default: 4000).

	MODULE 		: Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	SHELLCMD	: Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test")

```

### Parameters from command line arguments

For any parameters not defined in cmd. line arguments, default value will be used.

'bds $BDS_SCRIPT -h' shows help for all command line arguments.

```
	-c <string>       : Configuration file path (if not specified, define parameters in command line argument).

	-prefix <string>  : Prefix for all outputs.
	-o <string>       : Output directory. (default: out)
	-tmp <string>     : Temporary directory for intermediate files. (default: tmp).

	-wt <int>         : Default walltime in seconds for all cluster jobs (default: 36000).
	-nth <int>        : Default number of threads for all cluster jobs (default: 1).
	-mem <int>        : Default max. memory in MB for all cluster jobs (default: 4000).

	-mod <string>     : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	-shcmd <string>   : Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test")
```


### How to specify softwares for the pipeline?

With a configuration file, you can specify softwares with MODULE* and SHELLCMD*. Any suffix is allowed like MODULE_BIO, MODULE_LANG, SHELLCMD_TOOL, SHELLCMD_TEST.

For example, to use bwa 0.7.7 and bedtools 2.x.x and samtools 1.2

```
MODULE_BWA= bwa/0.7.7 
MODULE_ANY_SUFFIX_SHOULD_BE_OKAY= bedtools/2.x.x; samtools/1.2
```

Additional environment variables can be defined with SHELLCMD*

```
# There is a bug in BDS code.
# YOU NEED TO COVER ANY ENVIRONMENT VARS WITH CURLY BRACKETS ${} !!
SHELLCMD_ETC= export PATH="${PATH}:/usr/bin/example"
SHELLCMD_TEST= TEST_PATH="/usr/lib/example"; export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${TEST_PATH}"
SHELLCMD_YOU_CAN_DEFINE_ANY_CUSTOM_VAR= CUST_VAR=200
```

With cmd. line arguments, use -mod [] and -shcmd []. Use semicolon as a delimiter.

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
