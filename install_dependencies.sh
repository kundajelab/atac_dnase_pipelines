#!/bin/bash

############ install in conda env.
conda config --add channels r
conda config --add channels bioconda
conda config --add channels astro

conda create -n bds_atac --file requirements.txt -y
conda create -n bds_atac_py3 --file requirements_py3.txt -y


############ install additional packages
source activate bds_atac

CONDA_ACTIVATE=$(which activate)
CONDA_BIN=$(dirname $(which activate))
CONDA_EXTRA="$CONDA_BIN/extras"

mkdir -p $CONDA_EXTRA

function add_to_activate {
  for i in "${CONTENTS[@]}"; do
    if [ $(grep "$i" "$CONDA_ACTIVATE" | wc -l ) == 0 ]; then
      echo $i >> "$CONDA_ACTIVATE"
    fi
  done
}

### BDS
CONTENTS=("export PATH=\$PATH:\$HOME/.bds")
add_to_activate

### PICARDROOT
CONTENTS=("export PICARDROOT=$CONDA_BIN")
add_to_activate

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
cd $CONDA_BIN
CONTENTS=("export PATH=\$PATH:$CONDA_EXTRA/phantompeakqualtools")
add_to_activate

##### preseq ==2.0.2 
cd $CONDA_EXTRA
rm -rf preseq
git clone https://github.com/smithlabcode/preseq --recursive
cd preseq
git checkout tags/v2.0.2
make all
cd $CONDA_BIN
ln -s $CONDA_EXTRA/preseq/preseq preseq

### texlive ==2013
cd $CONDA_EXTRA
mkdir -p texlive
wget http://mirror.ctan.org/systems/texlive/tlnet/install-tl-unx.tar.gz -N
tar zxvf install-tl-unx.tar.gz
rm -f install-tl-unx.tar.gz
cd install-tl*
echo "selected_scheme scheme-basic" > texlive_profile.txt
echo "TEXDIR $CONDA_EXTRA/texlive/2015" >> texlive_profile.txt
echo "TEXMFCONFIG ~/.texlive2015/texmf-config" >> texlive_profile.txt
echo "TEXMFHOME ~/texmf" >> texlive_profile.txt
echo "TEXMFLOCAL $CONDA_EXTRA/texlive/texmf-local" >> texlive_profile.txt
echo "TEXMFSYSCONFIG $CONDA_EXTRA/texlive/2015/texmf-config" >> texlive_profile.txt
echo "TEXMFSYSVAR $CONDA_EXTRA/texlive/2015/texmf-var" >> texlive_profile.txt
echo "TEXMFVAR ~/.texlive2015/texmf-var" >> texlive_profile.txt
echo "binary_x86_64-linux 1" >> texlive_profile.txt
echo "collection-basic 1" >> texlive_profile.txt
echo "collection-latex 1" >> texlive_profile.txt
echo "in_place 0" >> texlive_profile.txt
echo "option_adjustrepo 1" >> texlive_profile.txt
echo "option_autobackup 1" >> texlive_profile.txt
echo "option_backupdir tlpkg/backups" >> texlive_profile.txt
echo "option_desktop_integration 1" >> texlive_profile.txt
echo "option_doc 1" >> texlive_profile.txt
echo "option_file_assocs 1" >> texlive_profile.txt
echo "option_fmt 1" >> texlive_profile.txt
echo "option_letter 0" >> texlive_profile.txt
echo "option_menu_integration 1" >> texlive_profile.txt
echo "option_path 0" >> texlive_profile.txt
echo "option_post_code 1" >> texlive_profile.txt
echo "option_src 1" >> texlive_profile.txt
echo "option_sys_bin /usr/local/bin" >> texlive_profile.txt
echo "option_sys_info /usr/local/info" >> texlive_profile.txt
echo "option_sys_man /usr/local/man" >> texlive_profile.txt
echo "option_w32_multi_user 1" >> texlive_profile.txt
echo "option_write18_restricted 1" >> texlive_profile.txt
echo "portable 0" >> texlive_profile.txt
./install-tl -profile texlive_profile.txt
cd $CONDA_BIN
PDFLATEX=$(find . -name 'pdflatex' | head -n 1)
ln -s $CONDA_BIN/$PDFLATEX pdflatex



source deactivate
