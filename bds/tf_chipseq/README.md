TF ChIP-Seq Pipeline
===============================================

TF ChIP-Seq pipeline is based on https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#.
Please take a look at ../README.md first.

### Parameters from configuration file

```
	PREFIX         		: Prefix for all outputs.
	OUTPUT_DIR     		: Output directory. (default: out)
	TMP_DIR        		: Temporary directory for intermediate files. (default: tmp).

	WALLTIME       		: Default walltime in seconds for all cluster jobs (default: 36000).
	NTHREADS       		: Default number of threads for all cluster jobs (default: 1).
	MEMORY         		: Default max. memory in MB for all cluster jobs (default: 4000).

	MODULE          	: Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	SHELLCMD        	: Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test")
        ADDPATH                 : Paths to be added to env. var. PATH separated by ; or :. (a quicker way to add PATH)


	QC_ONLY                 : Set it true to test-run and stop before peak calling, false: keep going through IDR (default: false).
	NUM_REP                 : # of replicates, set it only for qc = true. (default: 1).

	NTHREADS_BWA            : Number of threads for bwa aln (default: 4).
	BWA_INDEX_NAME          : Path for bwa index.
	BWA_ALN_PARAM           : Parameters for bwa align (default: "-q 5 -l 32 -k 2" ).

	MAPQ_THRESH             : MAPQ_THRESH (default: 30).
	NREADS                  : Parameter for NREADS (default. 15000000).
	NTHREADS_RUN_SPP        : Number of threads for run_spp.R (default: 4).

	CREATE_WIG              : Set it true to create wig (default: false).
	CREATE_BEDGRAPH         : Set it true to create bedgraph (default: false).
	CONVERT_TO_BIGWIG       : Set it true to convert bedgraph to bigwig signal track (default: false).

	UMAP_DIR                : Path for umap (for hg19, path for globalmap_k20tok54).
	SEQ_DIR                 : Dir. for sequence files (for hg19, dir where chr*.fa exist).
	CHROM_SIZES             : Path for chrom.sizes file for your sequence files.

	INPUT_FASTQ_REP1        : Path for input fastq for replicate 1 (single ended).
	INPUT_FASTQ_REP2        : Path for input fastq for replicate 2 (single ended).

	INPUT_FASTQ_REP1_PE1    : Path for input fastq for replicate 1 pair 1 (paired-end).
	INPUT_FASTQ_REP1_PE2    : Path for input fastq for replicate 1 pair 2 (paired-end).
	INPUT_FASTQ_REP2_PE1    : Path for input fastq for replicate 2 pair 1 (paired-end).
	INPUT_FASTQ_REP2_PE2    : Path for input fastq for replicate 2 pair 2 (paired-end).

	INPUT_FASTQ_CTL_REP1    : Path for control fastq for replicate 1 (single ended).
	INPUT_FASTQ_CTL_REP2    : Path for control fastq for replicate 2 (single ended).
	
	INPUT_FASTQ_CTL_REP1_PE1: Path for control fastq for replicate 1 pair 1 (paired-end).
	INPUT_FASTQ_CTL_REP1_PE2: Path for control fastq for replicate 1 pair 2 (paired-end).
	INPUT_FASTQ_CTL_REP2_PE1: Path for control fastq for replicate 2 pair 1 (paired-end).
	INPUT_FASTQ_CTL_REP2_PE2: Path for control fastq for replicate 2 pair 2 (paired-end).

	NPEAK                   : Parameter for -npeak in phantompeakqual tool run_spp.R (default: 300000).
	DUPE_REMOVED            : Set it true if dupes are removed when aligning (default: true).
	IDR_THRESH              : IDR thresh (default: 0.02).

```


### Parameters from command line arguments

```
	-c <string>             : Configuration file path (if not specified, define parameters in command line argument).

	-prefix <string>        : Prefix for all outputs.
	-o <string>             : Output directory. (default: out)
	-tmp <string>           : Temporary directory for intermediate files. (default: tmp).

	-wt <int>               : Default walltime in seconds for all cluster jobs (default: 36000).
	-nth <int>              : Default number of threads for all cluster jobs (default: 1).
	-mem <int>              : Default max. memory in MB for all cluster jobs (default: 4000).

	-mod <string>           : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	-shcmd <string>         : Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test").
        -addpath <string>       : Paths to be added to env. var. PATH separated by ; or :. (a quicker way to add PATH)


	-qc                     : Set it true to test-run and stop before peak calling, false: keep going through IDR (default: false).
	-num_rep <int>          : # of replicates, set it only for qc = true. (default: 1).

	-nth_bwa <int>          : Number of threads for bwa aln (default: 4).
	-idx_bwa <string>       : Path for bwa index.
	-param_bwa <string>     : Parameters for bwa align (default: "-q 5 -l 32 -k 2" ).

	-mapq_thresh <int>      : MAPQ_THRESH (default: 30).
	-nreads <int>           : Parameter for NREADS (default. 15000000).
	-nth_spp <int>          : Number of threads for run_spp.R (default: 4).

	-wig                    : Set it true to create wig (default: false).
	-bedgraph               : Set it true to create bedgraph (default: false).
	-bigwig                 : Set it true to convert bedgraph to bigwig signal track (default: false).

	-umap <string>          : Path for umap (for hg19, path for globalmap_k20tok54).
	-seq <string>           : Dir. for sequence files (for hg19, dir where chr*.fa exist).
	-chrsz <string>         : Path for chrom.sizes file for your sequence files.

	-fastq1 <string>        : Path for input fastq for replicate 1 (single ended).
	-fastq2 <string>        : Path for input fastq for replicate 2 (single ended).

	-ctl_fastq1 <string>    : Path for control fastq for replicate 1 (single ended).
	-ctl_fastq2 <string>    : Path for control fastq for replicate 2 (single ended).

	-fastq1_1 <string>      : Path for input fastq for replicate 1 pair 1 (paired-end).
	-fastq1_2 <string>      : Path for input fastq for replicate 1 pair 2 (paired-end).
	-fastq2_1 <string>      : Path for input fastq for replicate 2 pair 1 (paired-end).
	-fastq2_2 <string>      : Path for input fastq for replicate 2 pair 2 (paired-end).

	-ctl_fastq1_1 <string>  : Path for control fastq for replicate 1 pair 1 (paired-end).
	-ctl_fastq1_2 <string>  : Path for control fastq for replicate 1 pair 2 (paired-end).
	-ctl_fastq2_1 <string>  : Path for control fastq for replicate 2 pair 1 (paired-end).
	-ctl_fastq2_2 <string>  : Path for control fastq for replicate 2 pair 2 (paired-end).

	-npeak <int>            : Parameter for -npeak in phantompeakqual tool run_spp.R (default: 300000).
	-dup_rm                 : Set it true if dupes are removed when aligning (default: true).
	-idr_thresh <string>    : IDR thresh (default: 0.02).

```

### IGNORE STEPS BELOW IF YOU ARE WORKING ON OUR LAB CLUSTER!



### Local installation instruction for R-2.15.1 (for spp)

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

### Local installation for run_spp.R: Anshul's phantompeakqualtool (for spp)

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

### Local installation instruction for Python3 and packages (for Nathan's IDR)
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
Add the following lines to your $HOME/.bashrc and to your configuration file.
```
export PATH="${PATH}:${HOME}/usr/local/bin"; export PYTHONPATH="${HOME}/local/lib/python3.4/site-packages:${PYTHONPATH}"
```

### Local installation instruction for IDR

```
cd $HOME
git clone --recursive https://github.com/nboley/idr.git
cd idr
$HOME/usr/local/bin/python3.4 setup.py install --prefix=$HOME/local/
```
Add the following lines to your $HOME/.bashrc and to your configuration file.
```
export PATH="${HOME}/local/bin"
```

### Local installation instruction for Wiggler (for generating signal tracks from alignments)

<a href="https://code.google.com/p/align2rawsignal/">https://code.google.com/p/align2rawsignal/</a>

```
cd $HOME
wget https://align2rawsignal.googlecode.com/files/align2rawsignal.2.0.tgz
tar zxvf align2rawsignal.2.0.tgz
```
Add the following lines to your $HOME/.bashrc and to your configuration file.
```
export PATH="${HOME}/align2rawsignal/bin"
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
