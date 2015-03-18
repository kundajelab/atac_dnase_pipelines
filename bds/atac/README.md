ATAC Seq Pipeline
===

ATAC Seq pipeline is based on https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#.

Dependencies: BigDataScript (BDS)


### Configuration file

Modify $CONF_FILE (by default: conf_atac.txt) to have your own settings.

```
* PREFIX 			: Prefix for all output files
* OUTPUT_DIR 		: Output directory (both relative and absolute paths work)
* TMP_DIR 			: Temporary folder for intermediate files during bwa alignment
* USE_BGZIP			: Index BED type files (for visualization in a genome browser). Make sure bgzip and tabix installed add their path to MODULE_*.

* WALLTIME 			: default walltime for all jobs (in seconds)
* NTHREADS 			: default # of threads for all jobs
* MEMORY			: default max. memory for all jobs (in bytes)

* TRIM_ADAPTERS 	: Path for trimAdapters.py

* BOWTIE_IDX		: Path (prefix) of bowtie2 index files
* NTHREADS_BWT2		: # of threads for bwt2
* MEMORY_BWT2		: Max. memory limit for bwt2
* STACK_BWT2		: Stack size for bwt2

* MAPQ_THRESH		: MAPQ_THRESH
* JVM_OPTS			: Java VM additional options (eg. -Xmx4G)

* ADJUST_BED_TN5	: Path for adjustBedTN5.sh

* genomeSize  		: hs by default
* chrSize 	 		: Location of chrom.sizes file for your .fa

* MODULE_* 			: Freely name suffix and specify RHS, then BDS will run "module add RHS"
* EXPORT_* 			: Freely name suffix and specify RHS, then BDS will add env. variable to bash shell
```

### How to specify bioinformatics software version?

You can specify software versions with MODULE_* and EXPORT_* in the above configuration file. You can freely name suffix of MODULE_ or EXPORT_ then BDS will read all keys starting with MODULE_ or EXPORT_. Use ; for new line.

For example, to use bwa 0.7.7 and bedtools 2.x.x and samtools 1.2

```
MODULE_BWA= bwa/0.7.7 
MODULE_ANY_SUFFIX_SHOULD_BE_OKAY= bedtools/2.x.x; samtools/1.2
```

The above lines in the configuration file will execute the following:

```
module add bwa/0.7.7; module add bedtools/2.x.x; module add samtools/1.2

```

To use a software unmoduled, For example "sw_example" is on /usr/bin/example/, then add software path to environment variable PATH.

```
EXPORT_ETC= export PATH="${PATH}:/usr/bin/example"
EXPORT_TEST= TEST_PATH="/usr/lib/example"; export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${TEST_PATH}"
EXPORT_YOU_CAN_DEFINE_ANY_CUSTOM_VAR= CUST_VAR=200
```

The above lines in the configuration file will execute the following:

```
export PATH="${PATH}:/usr/bin/example"; TEST_PATH="/usr/lib/example"; export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${TEST_PATH}"; CUST_VAR=200
```

IMPORTANT!! There is a bug in BDS code. YOU NEED TO WRAP ANY ENVIRONMENT VARS IN RHS WITH CURVED BRACKETS ${} !!!!!


### Usage 

Firstly, copy bds cluster settings (bds.config.scg3) in $PIPELINE_DIR (path for atac.bds) to BDS directory. 
```
cp $PIPELINE_DIR/bds.config.scg3 $HOME/.bds/bds.config
vi $HOME/.bds/bds.config
```

Modify the following line in $HOME/.bds/bds.config. Due to BDS bug, you cannot use $HOME. Instead specify absolute path to your BDS directory (/home/your_account_name/.bds by default).

```
clusterRunAdditionalArgs = -v PATH=/home/leepc12/.bds
```

You can also specify your own cluster configuration file but I already provided one for SCG3.
```
bds -c your_own_bds.config any_script.bds
```

It is recommend to separate your working directory $WORK_DIR from $PIPELINE_DIR because BDS produces a lot of intermediate dir/files on $WORK_DIR. 

If you want to locally run jobs on your machine with both configuration files not specifed,

```
cd $WORK_DIR
bds $PIPELINE_DIR/atac.bds
```

You can specify your own configuration file (conf_atac.txt by default) with -c.

```
bds $PIPELINE_DIR/atac.bds -c $CONF_FILE
```

If you want to run your jobs on scg3 cluster, use -s sge.
```
bds -s sge $PIPELINE_DIR/atac.bds
```

If you run BDS script with -dryRun, it does not actually run the job, it compiles the script and list jobs to be executed.
```
bds -dryRun -s sge $PIPELINE_DIR/atac.bds
```

For advanced users, both $BDS_CONFIG and $CONF_FILE (chipseq configuration file) can be specified like the following. Make sure you don't change the order of arguments.

```
bds -c $BDS_CONFIG -dryRun -s sge $PIPELINE_DIR/atac.bds -c $CONF_FILE
```

For debugging, BDS does not remove temporary bash scripts with parameter -l, -d for debugging info.
```
bds -l -d -s sge $PIPELINE_DIR/atac.bds
```


### Installation insctruction for BigDataScript (BDS)

Official homepage for BDS is at <a href="http://pcingola.github.io/BigDataScript/download.html">http://pcingola.github.io/BigDataScript/download.html</a> and the git repo for BDS is at <a href="https://github.com/pcingola/BigDataScript.git">https://github.com/pcingola/BigDataScript.git</a>.

```
cd $HOME
wget -O bds_Linux.tgz https://github.com/pcingola/BigDataScript/releases/download/v0.999h/bds_Linux.tgz
tar zxvf bds_Linux.tgz
```

BDS is installed at $HOME/.bds/. You should add it to your PATH. Edit your $HOME/.bashrc and add the following:
```
export PATH=$PATH:$HOME/.bds/
```

### Local installation instruction for Wiggler

<a href="https://code.google.com/p/align2rawsignal/">https://code.google.com/p/align2rawsignal/</a>

```
cd $HOME
wget https://align2rawsignal.googlecode.com/files/align2rawsignal.2.0.tgz
tar zxvf align2rawsignal.2.0.tgz
```

Download MCR2010b.bin and install .

```
wget http://www.broadinstitute.org/~anshul/softwareRepo/MCR2010b.bin
./MCR2010b.bin -console

```
For human genomes, download UMAP from here <a href="http://www.broadinstitute.org/~anshul/projects/umap/">http://www.broadinstitute.org/~anshul/projects/umap/</a>

Add the following lines to your $HOME/.bashrc and to your configuration file.

```
MCRROOT=<MCR_ROOT>/v714
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
#LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64 # BUGGY!!!, make sure you comment this
MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}
XAPPLRESDIR=${MCRROOT}/X11/app-defaults
export LD_LIBRARY_PATH
export XAPPLRESDIR

```

```
export PATH=${PATH}:$HOME/align2rawsignal/bin

```

### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
