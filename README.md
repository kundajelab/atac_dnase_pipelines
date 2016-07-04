ATAC Seq Pipeline
===================================================


### Installation instruction

<a href="https://github.com/kundajelab/bds_atac/blob/master/INSTALL.md">General</a>

<a href="https://github.com/kundajelab/bds_atac/blob/master/INSTALL_Kundaje.md">Kundaje lab</a>

<a href="https://github.com/kundajelab/bds_atac/blob/master/INSTALL_SCG.md">SCG3/4</a>


### Usage

We recommend using BASH to run pipelines.

1) Define parameters in command line argument.

For general use, use the following command line: (for PE data set)
```
$ bds atac.bds -fastq1_1 [READ1] -fastq1_2 [READ2] -bwt2_idx [BOWTIE2_INDEX] \
-gensz [GENOMESIZE; hs for human, mm for mouse] -chrsz [CHR_SIZES_FILE] -adapter [ADAPTER_TO_BE_TRIMMED]
```

The pipeline can automatically detect adapters if you remove `-adapter` from your command line. (<b>You need to have an access to the git repo https://github.com/nboley/GGR_code</b>)

To use old adapter trimmers (`trim_galore` for SE and `trimAdapter.py` for PE):
```
-old_trimmer
```

If your fastqs are already trimmed, add the following to the command line to skip trimming stage.
```
-trimmed_fastq
```
<b>IMPORTANT!</b> If your data set is <b>SINGLE ENDED</b> add the following to the command line, otherwise the pipeline works for PE by default.
```
-se 
```
For ATAQC, define the following parameters. See help (`$ bds atac.bds`) for description of all parameters. Data files for running ataqc (hg19 and mm9) can be found here: http://mitra.stanford.edu/kundaje/dskim89/public/ataqc/. Even though you don't use a species file `-species_file`, you need to specify a species name for ATAQC. You will get an ATAQC report per replicate. ATAQC is avaible only when you start a pipeline with FASTQ inputs.
```
-species [hg19, mm9 or ...] -tss_enrich [] -ref_fa [] -blacklist [] -dnase [] -prom [] -enh [] -reg2map [] -roadmap_meta []
```
If you want to just align your data (no peak calling or further steps like IDR).
```
-align
```
If you don't want ATAQC, add the following to command line. 
```
-no_ataqc 
```
If you have just one replicate (PE), define fastqs with `-fastq[REP_ID]_[PAIR_ID]`.
```
-fastq1_1 [READ_PAIR1] -fastq1_2 [READ_PAIR2] \
```
For multiple replicates (PE), define fastqs with `-fastq[REP_ID]_[PAIR_ID]`. Add -fastq[]_[] for each replicate and pair to the command line:replicates.
```
-fastq1_1 [READ_REP1_PAIR1] -fastq1_2 [READ_REP1_PAIR2] -fastq2_1 [READ_REP2_PAIR1] -fastq2_2 [READ_REP2_PAIR2] ...
```
For multiple replicates (SE), define fastqs with `-fastq[REP_ID]`:
```
-se -fastq1 [READ_REP1] -fastq2 [READ_REP2] ...
```
You can also specify an adapter to be trimmed for each fastq. Example: `-adapter1_1 [ADAPTER1_1] -adapter1_2 [ADAPTER1_2] ...` for PE or `-adapter1 [ADAPTER1] -adapter2 [ADAPTER2] ...`. Define adapters just like how you defined your fastqs.

You can start from bam files. There are two kinds of bam files (raw or deduped) and you need to explicitly choose between raw bam (bam) and deduped one (filt_bam) with `-input [BAM_TYPE]`. Don't forget to add `-se` if they are not paired end (PE). For raw bams,
```
-bam1 [RAW_BAM_REP1] -bam2 [RWA_BAM_REP1] ...
```
For deduped (filtered) bams:
```
-filt_bam1 [FILT_BAM_REP1] -filt_bam2 [FILT_BAM_REP1] ...
```
For tagaligns (non-tn5-shifted):
```
-tag1 [TAGALIGN_REP1] -tag2 [TAGALIGN_REP2] ...
```
To subsample beds (tagaligns) add the following to the command line. This is different from subsampling for cross-corr. analysis. Peaks will be called with subsampled tagaligns.
```
-subsample [NO_READS_TO_SUBSAMPLE]
```
To change # of lines to subsample for cross-corr. analysis. This will not affect tasks downstream (peak calling and IDR).
```
-nreads [NO_READS_TO_SUBSAMPLE]
```
To disable pseudo replicate generation, add the following. By default, peak calling and IDR will be done for true replicates and pseudo replicates, but if you have `-true_rep` in the command line, you will also get IDR on true replicates only.
```
-true_rep
```
IDR analysis is included in the pipeline by default. For better IDR QC, define blacklist idr (for `hg19`, http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz). For other genomes, <a href="https://sites.google.com/site/anshulkundaje/projects/blacklists">https://sites.google.com/site/anshulkundaje/projects/blacklists</a>
```
-blacklist [BLACKLIST_BED]
```
If you don't want IDR, add the following:
```
-no_idr
```
To list all parameters and default values for them,
```
$ bds atac.bds
```


2) Define parameters in configuration file.
Key names in a configruation file are identical to parameter names on command line. 
```
$ bds atac.bds [CONF_FILE]

$ cat [CONF_FILE]
fastq1_1= [READ1]
fastq1_2= [READ2]
...
```


### Species file and Environment file

For Kundaje clusters and SCG3/4, do not define any species specific parameters like `-tss_enrich`, `-ref_fa`, `-blacklist`, `-dnase`, `-prom`, `-enh`, `-reg2map`, `-roadmap_meta`, `-bwt2_idx`, `-gensz`, `-chrsz`. They are already defined in `./species/kundaje.env` and `./species/scg3.env`. You just need to define the name of species `-species [SPECIES; hg19, mm9, ...]`.
```
$ bds atac.bds -species [hg19 or mm9] ...
```
See details <a href="https://github.com/kundajelab/TF_chipseq_pipeline/blob/master/README_PIPELINE.md" target=_blank>here</a>



### Parallelization and multi-threading (IMPORTANT!)

For completely serialized jobs, add the following. Individual jobs can still go multi-threaded.
```
-no_par
```
You can set up a limit for total # threads. Total # threads used by a pipeline will not exceed this limit. By default, it's 16 on SCG3/4 and 8 on Kundaje clusters and others.
```
-nth [MAX_TOTAL_NO_THREADS]
```
A pipeline automatically distributes `[MAX_TOTAL_NO_THREADS]` threads for jobs according to the corresponding input file sizes. For example of two fastqs (1GB and 2GB) with `-nth 6`, 2 and 4 threads are allocated for aligning 1GB and 2GB fastqs, respectively. The same policy applies to other multi-threaded tasks like deduping and peak calling.

However, all multi-threaded tasks (like bwa, bowtie2, spp and macs2) still have their own max. memory (`-mem_APPNAME [MEM_APP]`) and walltime (`-wt_APPNAME [WALLTIME_APP]`) settings. Max. memory is <b>NOT PER CPU</b>. On Kundaje cluster (with SGE flag activated `bds -s sge chipseq.bds ...`) or on SCG3/4, the actual shell command submitted by BDS for each task is like the following:
```
qsub -pe shm [NTH_ALLOCATED_FOR_APP] -h_vmem=[MEM_APP]/[NTH_ALLOCATED_FOR_APP] -h_rt=[WALLTIME_APP] ...
```
This ensures that total memory reserved for a cluster job equals to `[MEM_APP]`.



### How to efficiently manage multiple pipeline runs? (using UNIX screen)

`./utils/bds_scr` is a BASH script to create a detached screen for a BDS script and redirect stdout/stderr to a log file `[LOG_FILE_NAME]`. If a log file already exists, stdout/stderr will be appended to it. Monitor your pipeline with `tail -f [LOG_FILE_NAME]`. The only difference between `bds_scr` and `bds` is that you have `[SCR_NAME] [LOG_FILE_NAME]` between `bds` command and its parameters.
```
bds_scr [SCR_NAME] [LOG_FILE_NAME] atac.bds ...
```
You can skip `[LOG_FILE_NAME]` then a log file `[SCR_NAME].log` will be generated on the working directory.
```
bds_scr [SCR_NAME] atac.bds ...
```
You can also add any BDS parameters (like `-dryRun`, `-d` and `-s`). The following example is for running a pipeline on Sun Grid Engine.
```
bds_scr [SCR_NAME] [LOG_FILE_NAME] -s sge atac.bds ...
```
Once the pipeline run is done, the screen will be automatically closed. To kill a pipeline manually while it's running:
```
screen -X -S [SCR_NAME] quit
```
Or use `./utils/kill_scr`:
```
kill_scr [SCR_NAME]
```



### Temporary files on `$TMP` or `/tmp`

If you stop a BDS pipeline with `Ctrl+C` while calling peaks with `spp`. Temporary files generated by `Rscript` are not removed and they still exist on `$TMP` (or `/tmp` if not explicitly exported). You need to manually remove them.



### How to add more functions to pipeline source code?

<a href="https://github.com/kundajelab/TF_chipseq_pipeline/blob/master/README_CODE.md" target=_blank>https://github.com/kundajelab/TF_chipseq_pipeline/blob/master/README_CODE.md</a>



### Requirements 

For python2 (python 2.x >= 2.7) and R-3.x, <a href="https://github.com/kundajelab/bds_atac/blob/master/requirements.txt" target=_blank>here</a>
For python3, <a href="https://github.com/kundajelab/bds_atac/blob/master/requirements_py3.txt" target=_blank>here</a>



### Troubleshooting

See more troubleshootings <a href="https://github.com/kundajelab/TF_chipseq_pipeline/blob/master/README_PIPELINE.md" target=_blank>here</a>


1) samtools ncurses bug

Prepend a directory for `libncurses.so.5` to `LD_LIBRARY_PATH`. See `install_dependencies.sh` for solution.
```
samtools: symbol lookup error: /lib/x86_64-linux-gnu/libncurses.so.5: undefined symbol: _nc_putchar
```

2) Error: could not find environment: bds_atac

Unload any Anaconda Python modules. Remove locally installed Anaconda Python from your `$PATH`.

### Alternate Cloud-based Implementations

* The <a href="https://www.encodeproject.org/pipelines/" target=_blank>Encyclopedia of DNA Elements (ENCODE) Project</a> is in the process of adopting this pipeline for uniform processing of ENCODE ATAC-seq data. The <a href="https://github.com/ENCODE-DCC" target=_blank>official ENCODE implementation</a> by the ENCODE Data Coordination Center will be an exact mirror of our pipeline on <a href="https://www.dnanexus.com/" target=_blank>the DNAnexus cloud</a> (i.e. results will be exactly reproducible). Note that using this service requires a user to pay for cloud compute time.

* <a href="http://www.epinomics.co/" target=_blank>Epinomics</a> provides an independent, *free*, cloud-based pipeline implementation that adheres to the analysis protocol specifications of our pipeline. This implementation can be accessed at <a href="https://open.epigenomics.co/#/encode" target=_blank>https://open.epigenomics.co/#/encode</a>

### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Chuan Sheng Foo - PhD Student, Computer Science Dept., Stanford University
* Daniel Kim - MD/PhD Student, Biomedical Informatics Program, Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
