#!/bin/bash
# Stop on error
set -e

## conda environment name

ENV_NAME=bds_atac
ENV_NAME_PY3=bds_atac_py3

## install wiggler or not

INSTALL_WIGGLER_AND_MCR=0
INSTALL_GEM=0
INSTALL_PEAKSEQ=0

## install packages from official channels (bioconda and r)

conda create -n ${ENV_NAME} --file requirements.txt -y -c defaults -c bioconda -c r -c bcbio -c daler -c asmeurer
conda create -n ${ENV_NAME_PY3} --file requirements_py3.txt -y -c defaults -c bioconda -c r -c bcbio -c daler -c asmeurer

### bash function definition

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

## install useful tools for BigDataScript

mkdir -p $HOME/.bds
cp --remove-destination ./utils/bds_scr ./utils/bds_scr_5min ./utils/kill_scr bds.config $HOME/.bds/
cp --remove-destination -rf ./utils/clusterGeneric/ $HOME/.bds/

## install additional packages

source activate ${ENV_NAME}

conda uninstall graphviz -y # graphviz in bioconda has segmentation fault bug
conda install graphviz -c anaconda -y

conda install ucsc-bedgraphtobigwig -c bioconda -y
conda install ucsc-bedtobigbed -c bioconda -y

CONDA_BIN=$(dirname $(which activate))
CONDA_EXTRA="$CONDA_BIN/../extra"
CONDA_ACTIVATE_D="$CONDA_BIN/../etc/conda/activate.d"
CONDA_INIT="$CONDA_ACTIVATE_D/init.sh"
CONDA_LIB="$CONDA_BIN/../lib"
if [[ $(find $CONDA_LIB -name '*egg-info*' -not -perm -o+r | wc -l ) > 0 ]]; then
  find $CONDA_LIB -name '*egg-info*' -not -perm -o+r -exec dirname {} \; | xargs chmod o+r -R
fi

mkdir -p $CONDA_EXTRA $CONDA_ACTIVATE_D

### install Anshul's phantompeakqualtool
cd $CONDA_EXTRA
git clone https://github.com/kundajelab/phantompeakqualtools
chmod 755 -R phantompeakqualtools
CONTENTS=("export PATH=$CONDA_EXTRA/phantompeakqualtools:\$PATH")
add_to_activate

### disable locally installed python package lookup
CONTENTS=("export PYTHONNOUSERSITE=True")
add_to_activate

if [ ${INSTALL_WIGGLER_AND_MCR} == 1 ]; then
  conda install -y -c conda-forge bc
  ### install Wiggler (for generating signal tracks)
  cd $CONDA_EXTRA
  wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/align2rawsignal/align2rawsignal.2.0.tgz -N --no-check-certificate
  tar zxvf align2rawsignal.2.0.tgz
  rm -f align2rawsignal.2.0.tgz
  CONTENTS=("export PATH=\$PATH:$CONDA_EXTRA/align2rawsignal/bin")
  add_to_activate

  ### install MCR (560MB)
  cd $CONDA_EXTRA
  wget http://mitra.stanford.edu/kundaje/software/MCR2010b.bin -N --no-check-certificate
  #wget https://personal.broadinstitute.org/anshul/softwareRepo/MCR2010b.bin -N --no-check-certificate
  chmod 755 MCR2010b.bin
  echo '-P installLocation="'${CONDA_EXTRA}'/MATLAB_Compiler_Runtime"' > tmp.stdin
  ./MCR2010b.bin -silent -options "tmp.stdin"
  rm -f tmp.stdin
  rm -f MCR2010b.bin
  CONTENTS=(
  "MCRROOT=${CONDA_EXTRA}/MATLAB_Compiler_Runtime/v714" 
  "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRROOT}/runtime/glnxa64" 
  "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRROOT}/bin/glnxa64" 
  "MCRJRE=\${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64" 
  "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRJRE}/native_threads" 
  "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRJRE}/server" 
  "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRJRE}" 
  "XAPPLRESDIR=\${MCRROOT}/X11/app-defaults" 
  "export LD_LIBRARY_PATH" 
  "export XAPPLRESDIR")
  add_to_activate
fi

# install PeakSeq
if [ ${INSTALL_PEAKSEQ} == 1 ]; then
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

CONDA_BIN=$(dirname $(which activate))
CONDA_EXTRA="$CONDA_BIN/../extra"
mkdir -p $CONDA_EXTRA

### uninstall IDR 2.0.3 and install the latest one
conda uninstall idr -y
cd $CONDA_EXTRA
git clone https://github.com/kundajelab/idr
cd idr
python3 setup.py install
cd $CONDA_EXTRA
rm -rf idr

# install GEM
if [ ${INSTALL_GEM} == 1 ]; then
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
