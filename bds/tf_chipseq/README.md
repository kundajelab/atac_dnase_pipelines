TF ChIP-Seq Pipeline
===

TF ChIP-Seq pipeline is based on https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#.

Please take a look at ../README.md .


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
```
Local installation for run_spp.R: Anshul's phantompeakqualtool
```
cd $HOME
wget https://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz
tar zxvf ccQualityControl.v.1.1.tar.gz
```
Make sure that you change mod for *.R in SPP to allow linux which finds them
```
chmod 755 phantompeakqualtools/*
```
Install spp package
```
mkdir $HOME/RLib
$HOME/R/bin/R CMD INSTALL -l $HOME/RLib $HOME/phantompeakqualtools/spp_1.10.1.tar.gz
```
Add the following line to your $HOME/.bashrc and to your configuration file.
```
export PATH=${PATH}:${HOME}/phantompeakqualtools
export R_LIBS=${HOME}/RLib
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
