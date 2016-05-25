#!/bin/bash
# Stop on error
set -e

conda create -n bds_atac --file requirements.txt -y -c defaults -c bioconda -c r
conda create -n bds_atac_py3 --file requirements_py3.txt -y -c defaults -c bioconda -c r

############ install additional packages

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

source activate bds_atac

CONDA_BIN=$(dirname $(which activate))
CONDA_EXTRA="$CONDA_BIN/../extra"
CONDA_LIB="$CONDA_BIN/../lib"
CONDA_ACTIVATE_D="$CONDA_BIN/../etc/conda/activate.d"
CONDA_INIT="$CONDA_ACTIVATE_D/init.sh"
mkdir -p $CONDA_EXTRA $CONDA_ACTIVATE_D

### BDS
mkdir -p $HOME/.bds
cp --remove-destination ./utils/bds_scr ./utils/kill_scr bds.config $HOME/.bds/
CONTENTS=("export PATH=\$PATH:\$HOME/.bds")
add_to_activate

### PICARDROOT
CONTENTS=("export PICARDROOT=$CONDA_BIN")
add_to_activate
cd $CONDA_BIN
ln -s ../share/picard-1.126-3/picard.jar picard.jar

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

#### install run_spp.R (Anshul's phantompeakqualtool)
cd $CONDA_EXTRA
wget https://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz -N
tar zxvf ccQualityControl.v.1.1.tar.gz
rm -f ccQualityControl.v.1.1.tar.gz
chmod 755 -R phantompeakqualtools
CONTENTS=("export PATH=\$PATH:$CONDA_EXTRA/phantompeakqualtools")
add_to_activate

source deactivate


source activate bds_atac_py3

if [ $? != 0 ]; then
  echo Anaconda environment not found!
  exit
fi

CONDA_BIN=$(dirname $(which activate))
CONDA_EXTRA="$CONDA_BIN/../extra"
mkdir -p $CONDA_EXTRA

# uninstall IDR 2.0.3 and install the latest one
conda uninstall idr -y
cd $CONDA_EXTRA
git clone https://github.com/nboley/idr
cd idr
python3 setup.py install
cd $CONDA_EXTRA
rm -rf idr

source deactivate



echo === Installing dependencies successfully done. ===
