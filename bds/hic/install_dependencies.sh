#!/bin/bash

SOFTWARE=/srv/gsfs0/scratch/leepc12/software

cd $SOFTWARE
wget https://cran.r-project.org/src/base/R-3/R-3.2.0.tar.gz --no-check-certificate
tar zxvf R-3.2.0.tar.gz
rm -f R-3.2.0.tar.gz
cd R-3.2.0
#./configure --prefix=$HOME/R --with-readline=no --with-x=no --enable-R-static-lib
./configure --with-readline=no --with-x=no --enable-R-static-lib --enable-R-shlib
make
cd $SOFTWARE
echo > tmp.R
  echo 'install.packages("snow", repos="http://cran.us.r-project.org")' >> tmp.R
  echo 'install.packages("snowfall", repos="http://cran.us.r-project.org")' >> tmp.R
  echo 'install.packages("bitops", repos="http://cran.us.r-project.org")' >> tmp.R
  echo 'install.packages("caTools", repos="http://cran.us.r-project.org")' >> tmp.R
  echo 'install.packages("fdrtool", repos="http://cran.us.r-project.org")' >> tmp.R
#  echo 'source("http://bioconductor.org/biocLite.R")' >> tmp.R
#  echo 'biocLite("Rsamtools",suppressUpdates=TRUE)' >> tmp.R
$SOFTWARE/R-3.2.0/bin/Rscript tmp.R
rm -f tmp.R

export R_HOME=$SOFTWARE/R-3.2.0
export PATH=$R_HOME/bin:$PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$R_HOME/lib/

$SOFTWARE/python2.7/bin/pip2.7 install --install-option="--prefix=$SOFTWARE/python2.7" rpy2





sudo su

sudo apt-get install libbz2-dev
. /etc/profile.d/modules.sh
module add r/3.2.0
pip install rpy2




. /etc/profile.d/modules.sh
module add r/3.2.0
module load python_anaconda3/2.3.0

cp /lib/x86_64-linux-gnu/libreadline.so.6 /software/python_anaconda3/2.3.0/lib/

/software/python_anaconda3/2.3.0/bin/pip3 install rpy2
