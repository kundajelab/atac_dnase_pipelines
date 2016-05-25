### Installation instruction

BDS and all dependencies have already been installed on lab servers. Do not run `install_dependencies.sh` on these servers. Get the latest version of the pipeline.
```
$ git clone https://github.com/kundajelab/bds_atac --recursive
$ cd bds_atac
```
Replace BDS's default `bds.config` with a correct one:
```
$ mkdir -p $HOME/.bds
$ cp bds.config ./utils/bds_scr $HOME/.bds
```


### Conflicts with local Anaconda Python

If you have a local Anaconda Python in your `$PATH`, the pipeline will not work correctly because it cannot find conda environments `bds_atac` and `bds_atac_py3` on your local Anaconda Python. Remove your local Anaconda Python from your `$PATH` or unload Anaconda Python modules. Check if `conda` points to the correct global `Miniconda3` installed on `/software/miniconda3/bin/conda`.
```
which conda
```

If you want to keep using your local Anaconda Python, run the following to install dependencies on your local one (this is not recommended):
```
./install_dependencies.sh
```



### Genome data files

For Kundaje lab servers (mitra, nandi, durga, kali, vayu, amold, wotan and kadru), the pipeline automatically determines the type of servers and set shell environments and species database. Skip all genome-specific parameters (e.g. bwa index) and just specify species.
```
$ bds atac.bds ... -species [SPECIES; hg19, mm9, ... ]
```
