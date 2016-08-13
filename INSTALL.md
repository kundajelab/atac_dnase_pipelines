### Installation instruction (Java)

Install Java 8 (jdk >= 1.8 or jre >= 1.8) on your system. For Debian/Ubuntu (>14.10) based Linux: `$ sudo apt-get install git openjdk-8-jre`. For Fedora/Red-Hat based Linux: `$ sudo yum install git java-1.8.0-openjdk`. For Ubuntu 14.04:
```
$ sudo add-apt-repository ppa:webupd8team/java -y
$ sudo apt-get update
$ sudo apt-get install oracle-java8-installer
```
If you don't have super-user privileges on your system, locally install Java 8 and add it to `$PATH`.


### Installation instruction (conda)

Install Anaconda Python3 (or Miniconda3 4.0.5) on your system. Recent versions of conda (>4.0.10) is buggy in parallel activation and do not work correctly with BDS. If you already have your own conda, downgrade it to 4.0.5 and skip Miniconda3 installation: `conda install conda=4.0.5`.
Get the Miniconda3 installer (4.0.5) at <a href="https://repo.continuum.io/miniconda/Miniconda3-4.0.5-Linux-x86_64.sh" target=_blank>https://repo.continuum.io/miniconda/Miniconda3-4.0.5-Linux-x86_64.sh</a> and install it.
```
$ wget https://repo.continuum.io/miniconda/Miniconda3-4.0.5-Linux-x86_64.sh
$ bash Miniconda3-4.0.5-Linux-x86_64.sh
```
Choose `yes` for the final question. If you choose `no`, you need to manually add Miniconda3 to your `$HOME/.bashrc`.
```
Do you wish the installer to prepend the Miniconda3 install location
to PATH in your /your/home/.bashrc ? [yes|no]
[no] >>> yes
```
Remove any other Anaconda Python from your `$PATH`.
Check your loaded modules with `$ module list` and unload any Anaconda Python modules. Open a new terminal after installation.


### Installation instruction (BigDataScript)

Install BigDataScript v0.99999e on your system.
```
$ wget https://github.com/leepc12/BigDataScript/blob/master/distro/bds_Linux.tgz?raw=true -O bds_Linux.tgz
$ mv bds_Linux.tgz $HOME
$ cd $HOME && tar zxvf bds_Linux.tgz
```
Add `export PATH=$PATH:$HOME/.bds` to your bash initialization script (`$HOME/.bashrc` or `$HOME/.bash_profile`).
If Java memory occurs, add `export _JAVA_OPTIONS="-Xms256M -Xmx512M -XX:ParallelGCThreads=1"` too.


### Installation instruction (pipeline)

Get the latest version of the pipeline.
```
$ git clone https://github.com/kundajelab/bds_atac --recursive
```
Install software dependencies automatically. It will create two conda environments (bds_atac and bds_atac_py3) under your conda.
```
$ ./install_dependencies.sh
```
If you see the following error, then update your Anaconda with `conda update conda` and downgrade it to 4.0.5 `conda install conda==4.0.5`.
See <a href="https://github.com/kundajelab/TF_chipseq_pipeline/issues/8">here</a> for more detail.
```
Error: ERROR: placeholder '/root/miniconda3/envs/_build_placehold_placehold_placehold_placehold_placehold_p' too short in: glib-2.43.0-2
```
If you don't use `install_dependencies.sh`, manually replace BDS's default `bds.config` with a correct one:
```
$ cp bds.config ./utils/bds_scr $HOME/.bds
```
If `install_dependencies.sh` fails, run `./uninstall_dependencies.sh`, fix problems and then try `./install_dependencies.sh` again.


### Installation instruction (genome data)

Install genome data for a specific genome `[GENOME]`. Currently `hg19`, `mm9`, `hg38`(BETA) and `mm10`(BETA) are available. Specify a directory `[DATA_DIR]` to download/process all genome data. A species file generated on `[DATA_DIR]` will be automatically added to your `./default.env` so that the pipeline knows that you have installed genome data using `install_genome_data.sh`. If you want to install both genomes (`hg19` and `mm9`), make sure that you use the same directory `[DATA_DIR]` for them. Each genome data will be installed on `[DATA_DIR]/[GENOME]`. If you use other BDS pipelines, it is recommended to use the same directory `[DATA_DIR]` to save disk space.

`./install_genome_data.sh` can take longer than an hour for downloading data and building bowtie2 index.
```
$ ./install_genome_data.sh [GENOME] [DATA_DIR]
```


### How to install dependencies and share them with others

If you have super-user privileges on your system, it is recommended to install Miniconda3 on `/opt/miniconda3/` and share conda environment with others.
```
$ sudo su
$ ./install_dependencies.sh
$ chmod 755 -R /opt/miniconda3/  # if you get some annoying permission issues.
```
In order to make Miniconda3 accessible for all users, create an intialization script `/etc/profile.d/conda_init.sh`.
```
$ echo '#!/bin/bash' > /etc/profile.d/conda_init.sh
$ echo 'export PATH=$PATH:/opt/miniconda3/bin' >> /etc/profile.d/conda_init.sh
```


### How to install genome data and share them with others

If you have super-user privileges on your system, you can install genome data on `/data/bds_pipeline_genome_data` and share them with others.
```
$ sudo su
$ ./install_genome_data.sh [GENOME] /data/bds_pipeline_genome_data
```
You can find a species file `[SPECIES_FILE]` on `/data/bds_pipeline_genome_data` for each pipeline type. Then others can use the genome data by adding `-species_file [SPECIES_FILE_PATH]` to the pipeline command line. Or they need to add `species_file = [SPECIES_FILE_PATH]` to their `./default.env`.
