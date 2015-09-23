HiC Pipeline
===================================================

### Installation instruction

Please read this README first!
<a href="https://github.com/kundajelab/ENCODE_chipseq_pipeline/blob/master/README_PIPELINE.md">README_PIPELINE.md</a>


### How does HiC pipeline works?

There are two stages for HiC pipeline.

1) Mapping and sorting (hic_map.bds)
: Inputs are fastqs. Map, align and sort them to generate cleaned pairs.

2) HiC (hic.bds)
: Using cleaned pairs from the stage 1), perform HiC analysis.

Two stages share the same indices for librarie (lib) and replicates (rep).


### Using species file

On Kundaje lab cluster, add the following to the command line, then you don't have to define species specific parameters like bwa index, chromosome sizes file and umap.
```
-kundaje_lab -species hg19
```


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
-nth_bwa [# THREADS FOR BWA]
```

2) Define parameters in configuration file.
```
$ bds hic_map.bds [CONF_FILE]

$ cat [CONF_FILE]
libs = D3,D6,SC
reps = R2,R8
fastq_L[Lib_ID]_R[Rep_ID]_P[Pair_ID] = [FASTQ_FOR_LIB1_REP1_PAIR1]
...
bwa_idx = [BWA_INDEX]
nth_bwa = [# THREADS FOR BWA]
```


### Usage (HiC)

There are two methods (using RE_file or fixed windows). RE_file path must be defined for method 're'. DO NOT USE all-mappable umap (starting from 1bp) for HiC analysis. Detailed information about umap is found on <a href="https://github.com/kundajelab/ENCODE_chipseq_pipeline/blob/master/README_PIPELINE.md">README_PIPELINE.md</a>.

1) Using Root directory of mapped data
```
$ bds hic.bds \
-method [METHOD; re: using RE_file, fixed: fixed windows] \
-RE_file [RE_FILE_PATH] \
-libs D3,D6,SC \
-reps R2,R8 \
-root_mapped [ROOT_OF_MAPPED_DATA: OUT_DIR OF MAPPING STAGE] \
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
-method [METHOD; re: using RE_file, fixed: fixed windows] \
-RE_file [RE_FILE_PATH] \
-libs D3,D6,SC \
-reps R2,R8 \
-cln_pair_L1_R1 [CLEANED_PAIR_FOR_LIB1_REP1] \
...
-merge [MERGE_METHOD; 0:no_merge, 1:merge_replicates_only, 2:merge_libraries_and_replicates] \
-res [COMMA-SEPERATED RESOLUTIONS; Example: 100,200,1000] \
-blacklist [BLACKLIST_FILE] \
-umap [MAPPABILITY_DATA] \
-chrsz [CHR_SIZES_FILE]
```


### Contributors

* Jin wook Lee - PhD Student, Mechanical Engineering Dept., Stanford University
* Anshul Kundaje - Assistant Professor, Dept. of Genetics, Stanford University
