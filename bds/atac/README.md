ATAC Seq Pipeline
===

Please take a look at ../README.md .

### Configuration file

Modify $CONF_FILE (by default: conf_atac.txt) to have your own settings.

```
* PREFIX 			: Prefix for all output files
* OUTPUT_DIR 		: Output directory (both relative and absolute paths work)
* TMP_DIR 			: Temporary folder for intermediate files during bwa alignment
* USE_BGZIP			: Index BED type files (for visualization in a genome browser). Make sure bgzip and tabix installed add their path to MODULE_*.

* WALLTIME 			: default walltime for all jobs (in seconds)
* NTHREADS 			: default # of threads for all jobs
* MEMORY			: default max. memory for all jobs (in bytes)

* TRIM_ADAPTERS 	: Path for trimAdapters.py

* BOWTIE_IDX		: Path (prefix) of bowtie2 index files
* NTHREADS_BWT2		: # of threads for bwt2
* MEMORY_BWT2		: Max. memory limit for bwt2
* STACK_BWT2		: Stack size for bwt2

* MAPQ_THRESH		: MAPQ_THRESH
* JVM_OPTS			: Java VM additional options (eg. -Xmx4G)

* ADJUST_BED_TN5	: Path for adjustBedTN5.sh

* genomeSize  		: hs by default
* chrSize 	 		: Location of chrom.sizes file for your .fa

* MODULE_* 			: Freely name suffix and specify RHS, then BDS will run "module add RHS"
* EXPORT_* 			: Freely name suffix and specify RHS, then BDS will add env. variable to bash shell
```


### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
