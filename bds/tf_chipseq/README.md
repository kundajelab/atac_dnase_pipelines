TF ChIP-Seq Pipeline
===============================================

TF ChIP-Seq pipeline is based on https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#.
Please take a look at <a href="../README.md">../README.md</a> first.

### Parameters from configuration file

```
$bds tf_chipseq.bds [CONF_FILE]

$cat [CONF_FILE]

	PREFIX         		: Prefix for all outputs.
	OUTPUT_DIR     		: Output directory. (default: out)
	TMP_DIR        		: Temporary directory for intermediate files. (default: tmp).

	REF_GENOME		: Reference genome name for epigenome browser track generation (eg. hg19, hg18, mm9 or mm10)

	WALLTIME       		: Default walltime in seconds for all cluster jobs (default: 36000).
	NTHREADS       		: Default number of threads for all cluster jobs (default: 1).
	MEMORY         		: Default max. memory in MB for all cluster jobs (default: 4000).

	MODULE          	: Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	SHELLCMD        	: Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test")
	ADDPATH                 : Paths to be added to env. var. PATH separated by ; or :. (a quicker way to add PATH)


	QC_ONLY                 : Set it true to test-run and stop before peak calling, false: keep going through IDR (default: false).
	NUM_REP                 : # of replicates, set it only for qc = true. (default: 1).

	NTHREADS_BWA            : Number of threads for bwa aln (default: 4).
	BWA_INDEX_NAME          : Path for bwa index.
	BWA_ALN_PARAM           : Parameters for bwa align (default: "-q 5 -l 32 -k 2" ).

	MAPQ_THRESH             : MAPQ_THRESH (default: 30).
	NREADS                  : Parameter for NREADS (default. 15000000).

	CREATE_WIG              : Set it true to create wig (default: false).
	CREATE_BEDGRAPH         : Set it true to create bedgraph (default: false).
	CONVERT_TO_BIGWIG       : Set it true to convert bedgraph to bigwig signal track (default: false).

	UMAP_DIR                : Path for umap (for hg19, path for globalmap_k20tok54).
	SEQ_DIR                 : Dir. for sequence files (for hg19, dir where chr*.fa exist).
	CHROM_SIZES             : Path for chrom.sizes file for your sequence files.

	INPUT_FASTQ_REP1        : Path for input fastq for replicate 1 (single ended).
	INPUT_FASTQ_REP2        : Path for input fastq for replicate 2 (single ended).

	INPUT_FASTQ_REP1_PE1    : Path for input fastq for replicate 1 pair 1 (paired-end).
	INPUT_FASTQ_REP1_PE2    : Path for input fastq for replicate 1 pair 2 (paired-end).
	INPUT_FASTQ_REP2_PE1    : Path for input fastq for replicate 2 pair 1 (paired-end).
	INPUT_FASTQ_REP2_PE2    : Path for input fastq for replicate 2 pair 2 (paired-end).

	INPUT_FASTQ_CTL_REP1    : Path for control fastq for replicate 1 (single ended).
	INPUT_FASTQ_CTL_REP2    : Path for control fastq for replicate 2 (single ended).
	
	INPUT_FASTQ_CTL_REP1_PE1: Path for control fastq for replicate 1 pair 1 (paired-end).
	INPUT_FASTQ_CTL_REP1_PE2: Path for control fastq for replicate 1 pair 2 (paired-end).
	INPUT_FASTQ_CTL_REP2_PE1: Path for control fastq for replicate 2 pair 1 (paired-end).
	INPUT_FASTQ_CTL_REP2_PE2: Path for control fastq for replicate 2 pair 2 (paired-end).

	NTHREADS_RUN_SPP        : Number of threads for spp (run_spp.R) (default: 4).
	NPEAK                   : Parameter for -npeak in phantompeakqual tool run_spp.R (default: 300000).
	DUPE_REMOVED            : Set it true if dupes are removed when aligning (default: true).
	IDR_THRESH              : IDR thresh (default: 0.02).

```


### Parameters from command line arguments

```
$bds tf_chipseq.bds [OPTS_FOR_TF_CHIPSEQ]

#[OPTS_FOR_TF_CHIPSEQ] are like the following:

	-c <string>             : Configuration file path (if not specified, define parameters in command line argument).

	-prefix <string>        : Prefix for all outputs.
	-o <string>             : Output directory. (default: out)
	-tmp <string>           : Temporary directory for intermediate files. (default: tmp).

	-gen <string>           : Reference genome name for epigenome browser track generation (eg. hg19, hg18, mm9 or mm10)

	-wt <int>               : Default walltime in seconds for all cluster jobs (default: 36000).
	-nth <int>              : Default number of threads for all cluster jobs (default: 1).
	-mem <int>              : Default max. memory in MB for all cluster jobs (default: 4000).

	-mod <string>           : Modules separated by ; (example: "bowtie/2.2.4; bwa/0.7.7; picard-tools/1.92").
	-shcmd <string>         : Shell cmds separated by ;. Env. vars should be written as ${VAR} not as $VAR (example: "export PATH=${PATH}:/usr/test; VAR=test").
        -addpath <string>       : Paths to be added to env. var. PATH separated by ; or :. (a quicker way to add PATH)


	-qc                     : Set it true to test-run and stop before peak calling, false: keep going through IDR (default: false).
	-num_rep <int>          : # of replicates, set it only for qc = true. (default: 1).

	-nth_bwa <int>          : Number of threads for bwa aln (default: 4).
	-idx_bwa <string>       : Path for bwa index.
	-param_bwa <string>     : Parameters for bwa align (default: "-q 5 -l 32 -k 2" ).

	-mapq_thresh <int>      : MAPQ_THRESH (default: 30).
	-nreads <int>           : Parameter for NREADS (default. 15000000).

	-wig                    : Set it true to create wig (default: false).
	-bedgraph               : Set it true to create bedgraph (default: false).
	-bigwig                 : Set it true to convert bedgraph to bigwig signal track (default: false).

	-umap <string>          : Path for umap (for hg19, path for globalmap_k20tok54).
	-seq <string>           : Dir. for sequence files (for hg19, dir where chr*.fa exist).
	-chrsz <string>         : Path for chrom.sizes file for your sequence files.

	-fastq1 <string>        : Path for input fastq for replicate 1 (single ended).
	-fastq2 <string>        : Path for input fastq for replicate 2 (single ended).

	-ctl_fastq1 <string>    : Path for control fastq for replicate 1 (single ended).
	-ctl_fastq2 <string>    : Path for control fastq for replicate 2 (single ended).

	-fastq1_1 <string>      : Path for input fastq for replicate 1 pair 1 (paired-end).
	-fastq1_2 <string>      : Path for input fastq for replicate 1 pair 2 (paired-end).
	-fastq2_1 <string>      : Path for input fastq for replicate 2 pair 1 (paired-end).
	-fastq2_2 <string>      : Path for input fastq for replicate 2 pair 2 (paired-end).

	-ctl_fastq1_1 <string>  : Path for control fastq for replicate 1 pair 1 (paired-end).
	-ctl_fastq1_2 <string>  : Path for control fastq for replicate 1 pair 2 (paired-end).
	-ctl_fastq2_1 <string>  : Path for control fastq for replicate 2 pair 1 (paired-end).
	-ctl_fastq2_2 <string>  : Path for control fastq for replicate 2 pair 2 (paired-end).

	-nth_spp <int>          : Number of threads for spp (run_spp.R) (default: 4).
	-npeak <int>            : Parameter for -npeak in phantompeakqual tool run_spp.R (default: 300000).
	-dup_rm                 : Set it true if dupes are removed when aligning (default: true).
	-idr_thresh <string>    : IDR thresh (default: 0.02).

```

### Installation instruction for softwares for the tf chipseq pipeline

Look at <a href="./README_SW.md">./README_SW.md</a>


### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
