### Installation instruction

BDS and all dependencies have already been installed on lab servers. Do not run `install_dependencies.sh` on these servers. Get the latest version of the pipeline.
```
$ git clone https://github.com/kundajelab/bds_atac --recursive
$ cd bds_atac
```
Replace BDS's default `bds.config` with a correct one:
```
$ mkdir -p $HOME/.bds
$ cp bds.config bds_scr $HOME/.bds
```


### Genome data files

For Kundaje lab servers (mitra, nandi, durga, kali, vayu, amold, wotan and kadru), the pipeline automatically determines the type of servers and set shell environments and species database. Skip all genome-specific parameters (e.g. bwa index) and just specify species.
```
$ bds atac.bds ... -species [SPECIES; hg19, mm9, ... ]
```
