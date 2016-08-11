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
```
Error: ERROR: placeholder '/root/miniconda3/envs/_build_placehold_placehold_placehold_placehold_placehold_p' too short in: glib-2.43.0-2
```
If you don't use `install_dependencies.sh`, manually replace BDS's default `bds.config` with a correct one:
```
$ cp bds.config ./utils/bds_scr $HOME/.bds
```
If `install_dependencies.sh` fails, run `./uninstall_dependencies.sh`, fix problems and then try `./install_dependencies.sh` again.


### Genome data files

For SCG3/4 (carmack*, crick*, scg3*, scg4*) and Sherlock (sherlock*.stanford.edu) clusters, the pipeline automatically determines the type of servers and set shell environments and species database. Skip all genome-specific parameters (e.g. bwa index) and just specify species.
```
$ bds atac.bds ... -species [SPECIES; hg19, mm9, ... ]
```
