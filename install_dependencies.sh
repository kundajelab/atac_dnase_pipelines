#!/bin/bash
# Stop on error
set -e

############ install in conda env.
conda config --add channels r # r-2.15.3
conda config --add channels bioconda
conda config --add channels astro # r-2.15.3
conda config --add channels kalefranz # bug free readline
conda config --add channels asmeurer # libtool

conda create -n bds_atac --file requirements.txt -y
conda install -n bds_atac --file requirements_r2.txt -y --force # force install R-2.15.3

conda create -n bds_atac_py3 --file requirements_py3.txt -y


function add_to_activate {
  if [ ! -f $CONDA_INIT ]; then
    echo > $CONDA_INIT
  fi
  for i in "${CONTENTS[@]}"; do
    if [ $(grep "$i" "$CONDA_INIT" | wc -l ) == 0 ]; then
      echo $i >> "$CONDA_INIT"
    fi
  done
}

############ install additional packages
source activate bds_atac

if [ $? != 0 ]; then
  echo Anaconda environment not found!
  exit
fi

CONDA_BIN=$(dirname $(which activate))
CONDA_EXTRA="$CONDA_BIN/../extra"
CONDA_LIB="$CONDA_BIN/../lib"
CONDA_ACTIVATE_D="$CONDA_BIN/../etc/conda/activate.d"
CONDA_INIT="$CONDA_ACTIVATE_D/init.sh"
mkdir -p $CONDA_EXTRA $CONDA_ACTIVATE_D

## GRAPHVIZ
# graphviz on bioconda is buggy (GLIBC 2.14 error)
conda uninstall graphviz -y
conda install graphviz -c defaults --override-channels -y

### BDS
mkdir -p $HOME/.bds
cp bds_scr bds.config $HOME/.bds/
CONTENTS=("export PATH=\$PATH:\$HOME/.bds")
add_to_activate

### PICARDROOT
CONTENTS=("export PICARDROOT=$CONDA_BIN")
add_to_activate
cd $CONDA_BIN
ln -s ../share/picard-1.126-3/picard.jar picard.jar

### SAMTOOLS NCURSES BUG FIX
CONTENTS=("export LD_LIBRARY_PATH=$CONDA_LIB:\$LD_LIBRARY_PATH")
add_to_activate

##### preseq ==2.0.2 
cd $CONDA_EXTRA
rm -rf preseq
git clone https://github.com/smithlabcode/preseq --recursive
cd preseq
git checkout tags/v2.0.2
export CPLUS_INCLUDE_PATH=$CONDA_BIN/../include
export LIBRARY_PATH=$CONDA_BIN/../lib
make all
cd $CONDA_BIN
ln -s $CONDA_EXTRA/preseq/preseq preseq


#### install R 2.x.x packages
cd $CONDA_EXTRA
wget http://mitra.stanford.edu/kundaje/software/spp_1.13.tar.gz -N
echo > tmp.R
echo 'withCallingHandlers(install.packages("snow", repos="http://cran.us.r-project.org"),warning = function(w) quit(save = "no", status = 1, runLast = FALSE))' >> tmp.R
echo 'withCallingHandlers(install.packages("snowfall", repos="http://cran.us.r-project.org"),warning = function(w) quit(save = "no", status = 1, runLast = FALSE))' >> tmp.R
echo 'withCallingHandlers(install.packages("bitops", repos="http://cran.us.r-project.org"),warning = function(w) quit(save = "no", status = 1, runLast = FALSE))' >> tmp.R
echo 'withCallingHandlers(install.packages("caTools", repos="http://cran.us.r-project.org"),warning = function(w) quit(save = "no", status = 1, runLast = FALSE))' >> tmp.R
echo 'source("http://bioconductor.org/biocLite.R")' >> tmp.R
echo 'biocLite("Rsamtools",suppressUpdates=TRUE)' >> tmp.R
echo 'withCallingHandlers(install.packages("./spp_1.13.tar.gz"),warning = function(w) quit(save = "no", status = 1, runLast = FALSE))' >> tmp.R
export RHOME=$(R RHOME)
Rscript tmp.R
rm -f tmp.R spp_1.13.tar.gz
CONTENTS=("export RHOME=\$(R RHOME)")
add_to_activate

#### install run_spp.R (Anshul's phantompeakqualtool)
cd $CONDA_EXTRA
wget https://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz -N
tar zxvf ccQualityControl.v.1.1.tar.gz
rm -f ccQualityControl.v.1.1.tar.gz
chmod 755 -R phantompeakqualtools
CONTENTS=("export PATH=\$PATH:$CONDA_EXTRA/phantompeakqualtools")
add_to_activate

source deactivate



#source activate bds_atac_py3
#
#if [ $? != 0 ]; then
#  echo Anaconda environment not found!
#  exit
#fi
#
#CONDA_BIN=$(dirname $(which activate))
#CONDA_EXTRA="$CONDA_BIN/../extra"
#mkdir -p $CONDA_EXTRA
#
## uninstall IDR 2.0.2 and install beta
#conda uninstall idr -y
#cd $CONDA_EXTRA
#git clone https://github.com/nboley/idr
#cd idr
#python3 setup.py install
#cd $CONDA_EXTRA
#rm -rf idr
#
#source deactivate
