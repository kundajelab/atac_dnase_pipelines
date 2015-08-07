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
