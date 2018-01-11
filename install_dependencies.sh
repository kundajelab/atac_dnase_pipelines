#!/bin/bash
# Stop on error
set -e

## conda environment name

ENV_NAME=bds_atac
ENV_NAME_PY3=bds_atac_py3

INSTALL_GEM=0
INSTALL_PEAKSEQ=0

## install packages from official channels (bioconda and r)

conda create -n ${ENV_NAME} --file requirements.txt -y -c defaults -c bioconda -c r -c bcbio -c daler -c asmeurer
conda create -n ${ENV_NAME_PY3} --file requirements_py3.txt -y -c defaults -c bioconda -c r -c bcbio -c daler -c asmeurer

### bash function definition

function add_to_activate {
  if [[ ! -f $CONDA_INIT ]]; then
    echo > $CONDA_INIT
  fi
  for i in "${CONTENTS[@]}"; do
    if [[ $(grep "$i" "$CONDA_INIT" | wc -l ) == 0 ]]; then
      echo $i >> "$CONDA_INIT"
    fi
  done
}

## install useful tools for BigDataScript

mkdir -p $HOME/.bds
cp -f ./utils/bds_scr ./utils/bds_scr_5min ./utils/kill_scr bds.config $HOME/.bds/
cp -rf ./utils/clusterGeneric/ $HOME/.bds/

## install additional packages

source activate ${ENV_NAME}

conda uninstall graphviz -y # graphviz in bioconda has segmentation fault bug
conda install graphviz -c anaconda -y

conda install ucsc-bedgraphtobigwig -c bioconda -y
conda install ucsc-bedtobigbed -c bioconda -y

#CONDA_BIN=$(dirname $(which activate))/../envs/${ENV_NAME}/bin
#CONDA_BIN=$(dirname $(which activate))
CONDA_BIN=$(dirname $(which bedtools))
CONDA_EXTRA="$CONDA_BIN/../extra"
CONDA_ACTIVATE_D="$CONDA_BIN/../etc/conda/activate.d"
CONDA_INIT="$CONDA_ACTIVATE_D/init.sh"
CONDA_LIB="$CONDA_BIN/../lib"
if [[ $(find $CONDA_LIB -name '*egg-info*' -not -perm -o+r | wc -l ) > 0 ]]; then
  find $CONDA_LIB -name '*egg-info*' -not -perm -o+r -exec dirname {} \; | xargs chmod o+r -R
fi

mkdir -p $CONDA_EXTRA $CONDA_ACTIVATE_D

### install Anshul's phantompeakqualtool
echo $CONDA_EXTRA
cd $CONDA_EXTRA
git clone https://github.com/kundajelab/phantompeakqualtools
chmod 755 -R phantompeakqualtools
CONTENTS=("export PATH=$CONDA_EXTRA/phantompeakqualtools:\$PATH")
add_to_activate

### disable locally installed python package lookup
CONTENTS=("export PYTHONNOUSERSITE=True")
add_to_activate
#CONTENTS=("export PYTHONPATH=$CONDA_LIB/python2.7/site-packages:\$PYTHONPATH")
#add_to_activate

# install PeakSeq
if [[ ${INSTALL_PEAKSEQ} == 1 ]]; then
  cd $CONDA_EXTRA
  wget http://archive.gersteinlab.org/proj/PeakSeq/Scoring_ChIPSeq/Code/C/PeakSeq_1.31.zip -N --no-check-certificate
  unzip PeakSeq_1.31.zip
  rm -f PeakSeq_1.31.zip
  cd PeakSeq
  make
  chmod 755 bin/PeakSeq
  cd $CONDA_BIN
  ln -s $CONDA_EXTRA/PeakSeq/bin/PeakSeq
fi

source deactivate


source activate ${ENV_NAME_PY3}

#CONDA_BIN=$(dirname $(which activate))/../envs/${ENV_NAME_PY3}/bin
#CONDA_BIN=$(dirname $(which activate))
CONDA_BIN=$(dirname $(which bedtools))
CONDA_EXTRA="$CONDA_BIN/../extra"
CONDA_ACTIVATE_D="$CONDA_BIN/../etc/conda/activate.d"
CONDA_INIT="$CONDA_ACTIVATE_D/init.sh"
CONDA_LIB="$CONDA_BIN/../lib"
if [[ $(find $CONDA_LIB -name '*egg-info*' -not -perm -o+r | wc -l ) > 0 ]]; then
  find $CONDA_LIB -name '*egg-info*' -not -perm -o+r -exec dirname {} \; | xargs chmod o+r -R
fi

mkdir -p $CONDA_EXTRA $CONDA_ACTIVATE_D

### uninstall IDR 2.0.3 and install the latest one
conda uninstall idr -y
cd $CONDA_EXTRA
git clone --branch 2.0.4.1 https://github.com/kundajelab/idr
cd idr
python3 setup.py install
cd $CONDA_EXTRA
rm -rf idr

### disable locally installed python package lookup
CONTENTS=("export PYTHONNOUSERSITE=True")
add_to_activate
CONTENTS=("export PYTHONPATH=$CONDA_LIB/python3.5/site-packages:\$PYTHONPATH")
add_to_activate

# install GEM
if [[ ${INSTALL_GEM} == 1 ]]; then
  cd $CONDA_EXTRA
  wget http://groups.csail.mit.edu/cgs/gem/download/gem.v3.0.tar.gz -N --no-check-certificate
  tar zxvf gem.v3.0.tar.gz  
  rm -f gem.v3.0.tar.gz  
  cd gem
  chmod 755 gem.jar
  cd $CONDA_BIN
  ln -s $CONDA_EXTRA/gem/gem.jar
fi

source deactivate


echo == Installing dependencies has been successfully done. ==
