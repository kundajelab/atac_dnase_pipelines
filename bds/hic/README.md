HiC Pipeline
===================================================

### Installation instruction

Please read this README first!
<a href="../README.md">README.md</a>


### How does HiC pipeline works?

There are two stages for HiC pipeline.

1) Mapping and sorting (hic_map.bds)
: Inputs are fastqs. Map, align and sort them to generate cleaned pairs.

2) HiC (hic.bds)
: Using cleaned pairs from the stage 1), perform HiC analysis.

Two stages share the same indices for librarie (lib) and replicates (rep).



### Usage (mapping and sorting)

The following examples are for libraries (libs: D3, D6 and SC) and replicates (reps: R2 and R8).

1) Define parameters in command line argument. 
```
$ bds hic_map.bds \
-libs D3,D6,SC \
-reps R2,R8 \
-fastq_L[Lib_ID]_R[Rep_ID]_P[Pair_ID] [FASTQ_FOR_LIB1_REP1_PAIR1] \
...
-bwa_idx [BWA_INDEX]
```

For Kundaje lab cluster and SCG3, skip parameters (bwa_idx) and just specify species.
```
$ bds hic_map.bds ... -species [hg19 or mm9]
```

To change resource settings (# of processor, max memory and walltime) for bwa align, add the following to command line:
```
-nth_bwa_aln [NTHREADS_BWA_ALN] -mem_bwa_aln [MEMORY_BWA_ALN; e.g. 10G] -wt_bwa_aln [WALLTIME_BWA_ALN; e.g. 100h]
```


2) Define parameters in configuration file.

Key names in a configruation file are identical to parameter names on command line. 
```
$ bds hic_map.bds [CONF_FILE]

$ cat [CONF_FILE]
libs = D3,D6,SC
reps = R2,R8
...
```


### Usage (HiC)

There are two methods (using RE_file or fixed windows). RE_file path must be defined for method 're'. DO NOT USE all-mappable umap (starting from 1bp) for HiC analysis. Detailed information about umap is found on <a href="https://github.com/kundajelab/ENCODE_chipseq_pipeline/blob/master/README_PIPELINE.md">README_PIPELINE.md</a>.

1) Using Root directory of mapped data
```
$ bds hic.bds \
-root_mapped [ROOT_OF_MAPPED_DATA: OUT_DIR OF MAPPING STAGE] \
-method [METHOD; re: using RE_file, fixed: fixed windows] \
-RE_file [RE_FILE_PATH] \
-libs D3,D6,SC \
-reps R2,R8 \
-read_len [READ_LENGTH] \
...
-merge [MERGE_METHOD; 0:no_merge, 1:merge_replicates_only, 2:merge_libraries_and_replicates] \
-res [COMMA-SEPERATED RESOLUTIONS; Example: 100,200,1000] \
-blacklist [BLACKLIST_FILE] \
-umap [MAPPABILITY_DATA] \
-chrsz [CHR_SIZES_FILE]
```

2) Using paths for individual cleaned pairs
```
$ bds hic.bds \
-cln_pair_L1_R1 [CLEANED_PAIR_FOR_LIB1_REP1] \
...
-method [METHOD; re: using RE_file, fixed: fixed windows] \
...
```

For Kundaje lab cluster and SCG3, skip parameters (umap and chrsz) and just specify species.
```
$ bds hic.bds ... -species [hg19 or mm9]
```

For low resolution (<=50), you will need very long walltime (timeout) for hic step2 and step4. Add the following to the command line (for example: 200 hours for step2 and 150 hours for step4)
```
-wt_hic2 200h -wt_hic4 150h
```

### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
