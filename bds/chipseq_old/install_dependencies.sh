################ local installation ##################

# DO NOT CHANGE NAMING OF SOFTWARE DIRECTORY!
SOFTWARE=$HOME/software_bds
BASHRC=$HOME/add_to_bashrc

echo
echo "IMPORTANT!"
echo
echo "Please make sure that the following softwares are installed on your system."
echo "These command lines need super-user privilege to install them, or ask your system administrator."
echo "=============================================================================="
echo "sudo apt-get install build-essential zlib1g-dev ncurses-dev gfortran libboost-all-dev"
echo "sudo apt-get install openssl libssl-dev"
echo "sudo apt-get install libfreetype6-dev"
echo "sudo apt-get install liblapack-dev"
echo "sudo apt-get install git"
echo
echo "# Check if you already have java installed on your system"
echo "sudo apt-get install openjdk-7-jre"
echo
echo "# Execute the following lines if you don't have python 2.7 installed on your system"
echo "sudo add-apt-repository ppa:fkrull/deadsnakes"
echo "sudo apt-get update"
echo "sudo apt-get install python2.7"
echo "=============================================================================="
read -p "Press [Enter] key to continue..."


echo
echo "Add the following lines to your \$HOME/.bashrc or \$HOME/.bash_profile."
echo "Once installation is done, we recommend to add $BASHRC to one of them too."
echo "=============================================================================="
echo "# Java settings"
echo "export _JAVA_OPTIONS=\"-Xms256M -Xmx512M -XX:ParallelGCThreads=1\""
echo "export MAX_JAVA_MEM=\"8G\""
echo "export MALLOC_ARENA_MAX=4"
echo
echo "# BigDataScript settings"
echo "export PATH=\$PATH:\$HOME/.bds"
echo "=============================================================================="
read -p "Press [Enter] key to continue..."
echo


echo
echo "=============================================================================="
echo "Starting automatic installation for dependencies for ChIP-seq pipeline."
echo "Make sure you have enough disk space (at least 2GB) on your file system."
echo "All dependencies will be installed under $SOFTWARE."
echo "=============================================================================="
read -p "Press [Enter] key to continue..."
echo

mkdir -p $SOFTWARE
echo "# Add the following lines to your $HOME/.bashrc or $HOME/.bash_profile " > $BASHRC

# Local installation for BigDataScript (latest)
cd $HOME
wget https://github.com/pcingola/BigDataScript/blob/master/distro/bds_Linux.tgz?raw=true -O bds_Linux.tgz --no-check-certificate
tar zxvf bds_Linux.tgz
rm -f bds_Linux.tgz
echo "" >> $BASHRC
echo "# Path for BigDataScript (latest)" > $BASHRC
echo "export PATH=\$PATH:\$HOME/.bds/" >> $BASHRC

# Installation for chipseq pipeline (latest)
#cd $SOFTWARE
#git clone https://github.com/kundajelab/pipelines

# Local installation for bwa (0.7.10)
cd $SOFTWARE
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download -O bwa-0.7.10.tar.bz2
tar jxf bwa-0.7.10.tar.bz2
rm -f bwa-0.7.10.tar.bz2
cd bwa-0.7.10
make
echo "" >> $BASHRC
echo "# Path for bwa (0.7.10)" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/bwa-0.7.10" >> $BASHRC

# Local installation for samtools (0.1.19)
cd $SOFTWARE
wget http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download -O samtools-0.1.19.tar.bz2
tar jxf samtools-0.1.19.tar.bz2
rm -f samtools-0.1.19.tar.bz2
cd samtools-0.1.19
make
echo "" >> $BASHRC
echo "# Path for samtools (0.1.19)" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/samtools-0.1.19" >> $BASHRC

# Local installation for bedtools (2.19.1)
cd $SOFTWARE
wget http://pkgs.fedoraproject.org/repo/pkgs/BEDTools/bedtools-2.19.1.tar.gz/58de5377c3fb1bc1ab5a2620cf48f846/bedtools-2.19.1.tar.gz
tar zxvf bedtools-2.19.1.tar.gz
rm -f bedtools-2.19.1.tar.gz
cd bedtools2-2.19.1
make
echo "" >> $BASHRC
echo "# Path for bedtools (2.19.1)" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/bedtools2-2.19.1/bin" >> $BASHRC

# Local installation for UCSC tools
cd $SOFTWARE
mkdir ucsc_tools
cd ucsc_tools
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
chmod 755 *
echo "" >> $BASHRC
echo "# Path for UCSC tools" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/ucsc_tools" >> $BASHRC

# Local installation for PICARD tools (1.92)
cd $SOFTWARE
wget http://sourceforge.net/projects/picard/files/picard-tools/1.92/picard-tools-1.92.zip/download -O picard-tools-1.92.zip
unzip picard-tools-1.92.zip
rm -f picard-tools-1.92.zip
cd picard-tools-1.92
chmod 755 *
echo "" >> $BASHRC
echo "# Path for PICARd tools" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/picard-tools-1.92" >> $BASHRC
echo "export PICARDROOT=$SOFTWARE/picard-tools-1.92" >> $BASHRC

# Local installation for run_spp.R (Anshul's phantompeakqualtool)
cd $SOFTWARE
wget https://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz
tar zxvf ccQualityControl.v.1.1.tar.gz
rm -f ccQualityControl.v.1.1.tar.gz
chmod 755 phantompeakqualtools/*
echo "" >> $BASHRC
echo "# Path for run_spp.R (Anshul's phantompeakqualtool)" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/phantompeakqualtools" >> $BASHRC

# Local installation instruction for R (2.15.1) and relevant packages
cd $SOFTWARE
wget http://cran.r-project.org/src/base/R-2/R-2.15.1.tar.gz
tar zxvf R-2.15.1.tar.gz
rm -f R-2.15.1.tar.gz
cd R-2.15.1
#./configure --prefix=$HOME/R --with-readline=no --with-x=no --enable-R-static-lib
./configure --with-readline=no --with-x=no --enable-R-static-lib
make
cd $SOFTWARE
#wget http://compbio.med.harvard.edu/Supplements/ChIP-seq/spp_1.10.tar.gz
echo > tmp.R
  echo 'install.packages("snow", repos="http://cran.us.r-project.org")' >> tmp.R
  echo 'install.packages("snowfall", repos="http://cran.us.r-project.org")' >> tmp.R
  echo 'install.packages("bitops", repos="http://cran.us.r-project.org")' >> tmp.R
  echo 'install.packages("caTools", repos="http://cran.us.r-project.org")' >> tmp.R
#  echo 'install.packages("./spp_1.10.tar.gz")' >> tmp.R
  echo 'install.packages("./phantompeakqualtools/spp_1.10.1.tar.gz")' >> tmp.R
$SOFTWARE/R-2.15.1/bin/Rscript tmp.R
rm -f tmp.R
echo "" >> $BASHRC
echo "# Path for R (2.15.1)" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/R-2.15.1/bin" >> $BASHRC


# Local installation instruction for Python (3.4.3) and relevant packages (for Nathan Boley's IDR)
cd $SOFTWARE
wget https://www.python.org/ftp/python/3.4.3/Python-3.4.3.tgz
tar zxvf Python-3.4.3.tgz
rm -f Python-3.4.3.tgz
cd Python-3.4.3
./configure --prefix=$SOFTWARE/python3.4
make altinstall prefix=$SOFTWARE/python3.4 exec-prefix=$SOFTWARE/python3.4
wget http://cython.org/release/Cython-0.22.tar.gz
tar zxvf Cython-0.22.tar.gz
cd Cython-0.22
$SOFTWARE/python3.4/bin/python3.4 setup.py install --prefix=$SOFTWARE/python3.4
cd $SOFTWARE
$SOFTWARE/python3.4/bin/pip3.4 install --install-option="--prefix=$SOFTWARE/python3.4" numpy
$SOFTWARE/python3.4/bin/pip3.4 install --install-option="--prefix=$SOFTWARE/python3.4" matplotlib
echo "" >> $BASHRC
echo "# Path for run_spp.R (Anshul's phantompeakqualtool)" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/python3.4/bin" >> $BASHRC
echo "export PYTHONPATH=$SOFTWARE/python3.4/lib/python3.4/site-packages:\$PYTHONPATH">> $BASHRC

# Local installation instruction for Nathan Boley's IDR
cd $SOFTWARE
git clone --recursive https://github.com/nboley/idr.git
cd idr
$SOFTWARE/python3.4/bin/python3.4 setup.py install --prefix=$SOFTWARE/python3.4
echo "" >> $BASHRC
echo "# Path for IDR (written by Nathan Boley)" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/idr/bin" >> $BASHRC

# Local installation instruction for Anshul Kundaje's IDR
cd $SOFTWARE
wget https://sites.google.com/site/anshulkundaje/projects/idr/idrCode.tar.gz?attredirects=0 -O idrCode.tar.gz
tar zxvf idrCode.tar.gz
rm -f idrCode.tar.gz
echo "" >> $BASHRC
echo "# Path for IDR (written by Anshul Kundaje)" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/idrCode" >> $BASHRC

# Local installation instruction for Wiggler (for generating signal tracks)
cd $SOFTWARE
wget https://align2rawsignal.googlecode.com/files/align2rawsignal.2.0.tgz
tar zxvf align2rawsignal.2.0.tgz
rm -f align2rawsignal.2.0.tgz
echo "" >> $BASHRC
echo "# Path for Wiggler" >> $BASHRC
echo "export PATH=\$PATH:$SOFTWARE/align2rawsignal/bin" >> $BASHRC
wget http://www.broadinstitute.org/~anshul/softwareRepo/MCR2010b.bin
chmod 755 MCR2010b.bin
echo 1 > tmp.stdin
echo $SOFTWARE"/MATLAB_Compiler_Runtime" >> tmp.stdin
echo 1 >> tmp.stdin
echo 1 >> tmp.stdin
echo 3 >> tmp.stdin
./MCR2010b.bin -console < tmp.stdin
rm -f tmp.stdin
rm -f MCR2010b.bin
echo "" >> $BASHRC
echo "# Path for MCR2010b.bin" >> $BASHRC
echo "MCRROOT=$SOFTWARE/MATLAB_Compiler_Runtime/v714" >> $BASHRC
echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRROOT}/runtime/glnxa64" >> $BASHRC
echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRROOT}/bin/glnxa64" >> $BASHRC
echo "MCRJRE=\${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64" >> $BASHRC
echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRJRE}/native_threads" >> $BASHRC
echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRJRE}/server" >> $BASHRC
echo "LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:\${MCRJRE}" >> $BASHRC
echo "XAPPLRESDIR=\${MCRROOT}/X11/app-defaults" >> $BASHRC
echo "export LD_LIBRARY_PATH" >> $BASHRC
echo "export XAPPLRESDIR" >> $BASHRC

# WARNING
echo
echo
echo Done Installing all dependencies for ChIP-Seq pipeline
echo Created $BASHRC
echo
echo IMPORTANT!
echo
echo Please add $BASHRC to your $HOME/.bashrc or $HOME/.bash_profile


