TF ChIP-Seq Pipeline
===

TF ChIP-Seq pipeline is based on https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#.

Dependencies: BigDataScript


### Configuration file

Modify conf_tf_chipseq.txt to have your own settings.

* PRE_IDR 			: Set it true if you just want to compute QC score and stop before peak calling (default: false)
* NUM_REP			: # of replicates you want to test (default:2)

* PREFIX 			: Prefix for all output files
* OUTPUT_DIR 		: Output directory (absolute path works too)
* TMP_DIR 			: Temporary folder for intermediate files during bwa alignment

* USE_BGZIP			: Index BED type files (for visualization in a genome browser). Make sure bgzip and tabix installed.

* NTHREADS_BWA 		: # of threads for bwa aligment
* BWA_INDEX_NAME	: Prefix of bwa index files (eg. if you have bwa index files including hg19_Male.bwt, BWA_INDEX_NAME=hg19_Male)
* BWA_PARAM			: Parameters for bwa alignment (default: -q 5 -l 32 -k 2)

* MARKDUP 			: Dupe remover location
* MAPQ_THRESH		: MAPQ_THRESH

* NTHREADS_R		: # of threads for peak calling (spp)
* R_SCRIPT			: Specify if you have local installation of R (default: Rscript)
* RUN_SPP_DIR 		: Location of Anshul's pantompeakqualtools
* DUPE_REMOVED		: If true, use run_spp.nodups.R instead of run_spp.R
* NREADS 			: NREADS (default: 15000000)
* NPEAK 			: -npeak NPEAK in run_spp.R (default: 300000)

* PYTHON3 			: Specify if you have local installation of python3 (default: python3)
* IDR 				: Location of Nathan's IDR code
* IDR_THRESH	 	: Threshold for IDR (default=0.02)

* CREATE_WIG  		: Create wig file from .tagAlign.gz
* CREATE_BEDGRAPH 	: Create bedGraph file from .tagAlign.gz
* CONVERT_TO_BIGWIG : Convert bedGraph to bigwig
* CHROM_SIZES 		: Location of chrom.sizes file for your .fa
* UMAP_DIR 			: Location of umap (for hg19, globalmap_k20tok54)
* SEQ_DIR 			: Location of sequence .fa files (for hg19, chr?.fa)


### Usage 

 
Make sure you have the configuration file conf_tf_chipseq.txt on your WORKING directory $WORK_DIR (not on your script directory $PIPELINE_DIR). If bds.config is not specified with -c {bds.config}, $HOME/.bds/bds.config will be used by default. bds.config for scg3 (based on SGE) is included in the git repo. Separate your working directory ($WORK_DIR) from $PIPELINE_DIR cause BDS produces a lot of intermediate dir/files on $WORK_DIR.

If you want to locally run jobs on your machine,

```
cd $WORK_DIR
bds $PIPELINE_DIR/tf_chipseq.bds
```

On scg3 cluster, specify bds.config for scg3.
```
bds -c bds.config.scg3 -s sge $PIPELINE_DIR/tf_chipseq.bds 
```
Parameter -c {bds.config} doesn't work if it is put after the BDS script on the command line (bug in BDS). 

If you run BDS script with -dryRun, it does not actually run the job, it compiles the script and list jobs to be executed.

```
bds -c bds.config.scg3 -dryRun -s sge $PIPELINE_DIR/tf_chipseq.bds 
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
Add the following line in your $HOME/.bashrc .
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

Add the following lines in your $HOME/.bashrc .

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

Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
