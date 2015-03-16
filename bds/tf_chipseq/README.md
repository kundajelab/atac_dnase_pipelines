TF ChIP-Seq Pipeline
===

TF ChIP-Seq pipeline is based on https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#.

Dependencies: BigDataScript (BDS)


### Configuration file

Modify $CONF_FILE (by default: conf_tf_chipseq.txt) to have your own settings.

```
* PREFIX 			: Prefix for all output files
* OUTPUT_DIR 		: Output directory (both relative and absolute paths work)
* TMP_DIR 			: Temporary folder for intermediate files during bwa alignment
* USE_BGZIP			: Index BED type files (for visualization in a genome browser). Make sure bgzip and tabix installed add their path to MODULE_*.

* WALLTIME 			: default walltime for all jobs (in seconds)
* NTHREADS 			: default # of threads for all jobs
* MEMORY			: default max. memory for all jobs (in bytes)

* PRE_IDR 			: Set it true if you just want to compute QC score and stop before peak calling (default: false)
* NUM_REP			: # of replicates you want to test (default:2)

* NTHREADS_BWA 		: # of threads for bwa aligment
* BWA_INDEX_NAME	: Prefix of bwa index files (eg. if you have bwa index files including hg19_Male.bwt, BWA_INDEX_NAME=hg19_Male)
* BWA_PARAM			: Parameters for bwa alignment (default: -q 5 -l 32 -k 2)

* MARKDUP 			: Dupe marker path
* MAPQ_THRESH		: MAPQ_THRESH

* NTHREADS_R		: # of threads for peak calling (spp)
* DUPE_REMOVED		: If true, use run_spp.nodups.R instead of run_spp.R
* NREADS 			: NREADS (default: 15000000)
* NPEAK 			: -npeak NPEAK in run_spp.R (default: 300000)

* IDR_THRESH	 	: Threshold for IDR (default=0.02)

* CREATE_WIG  		: Create wig file from .tagAlign.gz
* CREATE_BEDGRAPH 	: Create bedGraph file from .tagAlign.gz
* CONVERT_TO_BIGWIG : Convert bedGraph to bigwig
* CHROM_SIZES 		: Location of chrom.sizes file for your .fa
* UMAP_DIR 			: Location of umap (for hg19, globalmap_k20tok54)
* SEQ_DIR 			: Location of sequence .fa files (for hg19, chr?.fa)

* MODULE_* 			: Freely name suffix and specify RHS, then BDS will run "module add RHS"
* EXPORT_* 			: Freely name suffix and specify RHS, then BDS will add env. variable to bash shell
```

### How to specify bioinformatics software version?

You can specify software versions with MODULE_* and EXPORT_* in the above configuration file. You can freely name suffix of MODULE_ or EXPORT_ then BDS will read all keys starting with MODULE_ or EXPORT_. Use ; for new line.

For example, to use bwa 0.7.7 and bedtools 2.x.x and bedtools2 1.x.x.

```
MODULE_BWA= bwa/0.7.7 
MODULE_BEDBED_ANY_SUFFIX_SHOULD_BE_OKAY= bedtools/2.x.x; bedtools2/1.x.x
```

The above lines in the configuration file will execute the following:

```
module add bwa/0.7.7; module add bedtools/2.x.x; module add bedtools2/1.x.x

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

```
EXPORT_BASH= export PATH="${PATH}:/bin:/usr/bin"
```

IMPORTANT!! Do not remove the above line in the conf. file. Common linux commands like rm will not work without it.


### Usage 


There are two configuration files for 1) pipeline ($CONF_FILE) and 2) sge cluster ($BDS_CONFIG). If 1) is not specified, you need to put conf_tf_chipseq.txt in your working directory ($WORK_DIR). If 2) is not specified, $HOME/.bds/bds.config will be used by default. It is recommend to separate $WORK_DIR from $PIPELINE_DIR because BDS produces a lot of intermediate dir/files on $WORK_DIR.

If you want to locally run jobs on your machine with both configuration files not specifed,

```
cd $WORK_DIR
bds $PIPELINE_DIR/tf_chipseq.bds
```

You can specify the location of configuration file 1).

```
bds $PIPELINE_DIR/tf_chipseq.bds -c $CONF_FILE
```

If you want to run your jobs on scg3 cluster, specify 2) cluster configuration file for scg3. $BDS_CONFIG for scg3 (bds.config.scg3) is provided in the git repo.
```
$BDS_CONFIG = bds.config.scg3
bds -c $BDS_CONFIG -s sge $PIPELINE_DIR/tf_chipseq.bds -c $CONF_FILE
```

Make sure you don't change the order of arguments. If you run BDS script with -dryRun, it does not actually run the job, it compiles the script and list jobs to be executed.

```
bds -c $BDS_CONFIG -dryRun -s sge $PIPELINE_DIR/tf_chipseq.bds -c $CONF_FILE
```

### BigDataScript (BDS)

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

### Local installation instruction for R-2.15.1

```
cd $HOME
wget http://cran.r-project.org/src/base/R-2/R-2.15.1.tar.gz
tar zxvf R-2.15.1.tar.gz
cd R-2.15.1
./configure --prefix=$HOME/R --with-readline=no --with-x=no --enable-R-static-lib
make
make install
cd ..
wget http://compbio.med.harvard.edu/Supplements/ChIP-seq/spp_1.10.tar.gz
$HOME/R/bin/R
	install.packages("snow")
	install.packages("snowfall")
	install.packages("bitops")
	install.packages("caTools")
q()
mkdir $HOME/RLib
$HOME/R/bin/R CMD INSTALL -l ~/RLib spp_1.10.tar.gz
```
Add the following line to your $HOME/.bashrc and to your configuration file.
```
export R_LIBS=$HOME/RLib
```

### Local installation instruction for Python3 and packages
```
cd $HOME

echo Python3
wget https://www.python.org/ftp/python/3.4.3/Python-3.4.3.tgz
tar zxvf Python-3.4.3.tgz
cd Python-3.4.3
./configure --prefix=$HOME/usr/local
make altinstall prefix=$HOME/usr/local exec-prefix=$HOME/usr/local

echo Cython 
wget http://cython.org/release/Cython-0.22.tar.gz
cd Cython-0.22
$HOME/usr/local/bin/python3.4 setup.py install --prefix=$HOME/local/

echo Packages
cd $HOME
$HOME/usr/local/bin/pip3.4 install --install-option="--prefix=$HOME/local" numpy
$HOME/usr/local/bin/pip3.4 install --install-option="--prefix=$HOME/local" matplotlib
```

### Local installation instruction for IDR

```
cd $HOME
git clone --recursive https://github.com/nboley/idr.git
cd idr
$HOME/usr/local/bin/python3.4 setup.py install --prefix=$HOME/local/
```

### Local installation instruction for Wiggler

<a href="https://code.google.com/p/align2rawsignal/">https://code.google.com/p/align2rawsignal/</a>

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

### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
