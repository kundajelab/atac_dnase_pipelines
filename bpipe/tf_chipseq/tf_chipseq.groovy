NTHREADS=4 // bwa (align_bwa), Rscript (xcor)
TMP="/srv/gsfs0/scratch/leepc12"
BWA_INDEX_NAME="/srv/gsfs0/projects/kundaje/commonRepository/indexes/bwa_indexes/encodeHg19Male/v0.7.10/encodeHg19Male_bwa-0.7.10.fa"
BWA_PARAM="-q 5 -l 32 -k 2"
MARKDUP="/srv/gs1/software/picard-tools/1.92/MarkDuplicates.jar"
MAPQ_THRESH=30
RSCRIPT="/srv/gs1/software/R/R-2.15.1/bin/Rscript"
NREADS=15000000
NREADS_PER_MILLION=15
NPEAK=300000
SPEAK=220

load "func.groovy"


prepare = { // temporary initial step due to bpipe bug. without this, the second square bracket does not work well
	outputs = inputs 
}

align_bwa = { // inputs : two fastq files (paired end)

	// check if paired end
        branch.PAIRED_END = inputs.split().join().count("PE2") > 0 // check if "PE2" exists in concatenation of input strings

	// declaration
	branch.OFPREFIX = input1.replaceAll(".PE1","").replaceAll(".PE2","").replaceAll(".fastq.gz","")

	doc "align_bwa: ${OFPREFIX}"

	branch.FASTQ_FILE_1=input1
	if ( PAIRED_END ) branch.FASTQ_FILE_2=input2

	branch.SAI_FILE_1="${OFPREFIX}_1.sai"
	if ( PAIRED_END ) branch.SAI_FILE_2="${OFPREFIX}_2.sai"

	branch.RAW_SAM_FILE="${OFPREFIX}.raw.sam.gz"

	branch.RAW_BAM_PREFIX="${OFPREFIX}.raw.srt"
	branch.RAW_BAM_FILE="${RAW_BAM_PREFIX}.bam" //# To be stored                                                                                                                          
	branch.BADCIGAR_FILE="${TMP}/badReads${OFPREFIX}.tmp"
	branch.RAW_BAM_FILE_MAPSTATS="${RAW_BAM_PREFIX}.flagstat.qc" //# QC File
	
	// run for outputs
	produce( RAW_BAM_FILE, RAW_BAM_FILE_MAPSTATS ) {

		if ( PAIRED_END ) {
			multi 	"bwa aln ${BWA_PARAM} -t ${NTHREADS} ${BWA_INDEX_NAME} ${FASTQ_FILE_1} > ${SAI_FILE_1}",
				"bwa aln ${BWA_PARAM} -t ${NTHREADS} ${BWA_INDEX_NAME} ${FASTQ_FILE_2} > ${SAI_FILE_2}"
			exec "bwa sampe ${BWA_INDEX_NAME} ${SAI_FILE_1} ${SAI_FILE_2} ${FASTQ_FILE_1} ${FASTQ_FILE_2} | gzip -c > ${RAW_SAM_FILE}"
			exec "rm ${SAI_FILE_1} ${SAI_FILE_2}"

			//# ==============================================================                                         
			//# Remove read pairs with bad CIGAR strings and sort by position                                                                
			//# ==============================================================
	
			//# Find bad CIGAR read names                                                                                                                                                  
			exec """
				zcat ${RAW_SAM_FILE} | awk 'BEGIN {FS="\t" ; OFS="\t"} ! /^@/ && \$6!="*" { cigar=\$6; gsub("[0-9]+D","",cigar); n = split(cigar,vals,"[A-Z]"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length(\$10) ; if (s!=seqlen) print \$1"t"; }' | sort | uniq > ${BADCIGAR_FILE}
			"""
			
			//# Remove bad CIGAR read pairs                                                                                                                                                
			/*
			exec """
				if [[ \$(cat ${BADCIGAR_FILE} | wc -l) -gt 0 ]]
				then
				    zcat ${RAW_SAM_FILE} | grep -v -F -f ${BADCIGAR_FILE} | samtools view -Su - | samtools sort - ${RAW_BAM_PREFIX}
				else
				    samtools view -Su ${RAW_SAM_FILE} | samtools sort - ${RAW_BAM_PREFIX}
				fi
			"""
			*/
			if ( get_line_cnt( file(BADCIGAR_FILE) ) > 0 )
			    exec "zcat ${RAW_SAM_FILE} | grep -v -F -f ${BADCIGAR_FILE} | samtools view -Su - | samtools sort - ${RAW_BAM_PREFIX}"
			else
			    exec "samtools view -Su ${RAW_SAM_FILE} | samtools sort - ${RAW_BAM_PREFIX}"
			
			exec "rm ${BADCIGAR_FILE} ${RAW_SAM_FILE}"
			exec "samtools flagstat ${RAW_BAM_FILE} > ${RAW_BAM_FILE_MAPSTATS}"
		}
		else {
			exec "bwa aln ${BWA_PARAM} -t ${NTHREADS} ${BWA_INDEX_NAME} ${FASTQ_FILE_1} > ${SAI_FILE_1}",
			exec "bwa samse ${BWA_INDEX_NAME} ${SAI_FILE_1} ${FASTQ_FILE_1} | samtools view -Su - | samtools sort - ${RAW_BAM_PREFIX}"
			exec "rm ${SAI_FILE_1}"
			exec "samtools flagstat ${RAW_BAM_FILE} > ${RAW_BAM_FILE_MAPSTATS}"
		}
	}
}

post_align_filt = {

	doc "post_align_filt: ${OFPREFIX}"

	branch.FILT_BAM_PREFIX="${OFPREFIX}.filt.srt"
	branch.FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"

	if ( PAIRED_END ) {
		branch.TMP_FILT_BAM_PREFIX="tmp.${FILT_BAM_PREFIX}.nmsrt"
		branch.TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"
	}
	branch.TMP_FILT_BAM_FILE2="${FILT_BAM_PREFIX}.dupmark.bam"
	branch.DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc"

	branch.FINAL_BAM_PREFIX="${OFPREFIX}.filt.srt.nodup"
	branch.FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" //# To be stored
	branch.FINAL_BAM_INDEX_FILE="${FINAL_BAM_PREFIX}.bai"
	branch.FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" //# QC file
	branch.FINAL_NMSRT_BAM_PREFIX="${OFPREFIX}.filt.nmsrt.nodup"
	branch.FINAL_NMSRT_BAM_FILE="${FINAL_NMSRT_BAM_PREFIX}.bam" //# To be stored

	branch.PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

	produce( FINAL_BAM_FILE, FINAL_BAM_INDEX_FILE, FINAL_BAM_FILE_MAPSTATS, PBC_FILE_QC ) {
		//# =============================
		//# Remove  unmapped, mate unmapped
		//# not primary alignment, reads failing platform
		//# Remove low MAPQ reads
		//# Only keep properly paired reads
		//# Obtain name sorted BAM file
		//# ==================
		if ( PAIRED_END ) {
			exec "samtools view -F 1804 -f 2 -q ${MAPQ_THRESH} -u ${RAW_BAM_FILE} | samtools sort -n - ${TMP_FILT_BAM_PREFIX} # Will produce name sorted BAM"
			
			//# Remove orphan reads (pair was removed)
			//# and read pairs mapping to different chromosomes
			//# Obtain position sorted BAM
			
			exec "samtools fixmate -r ${TMP_FILT_BAM_FILE} - | samtools view -F 1804 -f 2 -u - | samtools sort - ${FILT_BAM_PREFIX} # Will produce coordinate sorted BAM"
	
			exec "rm ${TMP_FILT_BAM_FILE}"		

			//# =============
			//# Mark duplicates
			//# =============
		
			exec "java -Xmx4G -jar ${MARKDUP} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE2} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
			exec "mv ${TMP_FILT_BAM_FILE2} ${FILT_BAM_FILE}"
		
			//# ============================
			//# Remove duplicates
			//# Index final position sorted BAM
			//# Create final name sorted BAM
			//# ============================

			exec "samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}"
			exec "samtools sort -n ${FINAL_BAM_FILE} ${FINAL_NMSRT_BAM_PREFIX}"
	
			//# Index Final BAM file
			exec "samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}"
			exec "samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}"

			//# =============================
			//# Compute library complexity
			//# =============================
			//# Sort by name
			//# convert to bedPE and obtain fragment coordinates
			//# sort by position and strand
			//# Obtain unique count statistics
		
			//# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
		
			exec """		
				samtools sort -no ${FILT_BAM_FILE} - | bamToBed -bedpe -i stdin | awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$4,\$6,\$9,\$10}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+\$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
			"""
			exec "rm ${FILT_BAM_FILE}"
		}
		else {
			exec "samtools view -F 1804 -q ${MAPQ_THRESH} -b ${RAW_BAM_FILE} > ${FILT_BAM_FILE}"
			//# ========================
			//# Mark duplicates
			//# ======================

			exec " java -Xmx4G -jar ${MARKDUP} INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE2} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
			exec "mv ${TMP_FILT_BAM_FILE2} ${FILT_BAM_FILE}"

			//# ============================
			//# Remove duplicates
			//# Index final position sorted BAM
			//# ============================

			exec "samtools view -F 1804 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}"

			//# Index Final BAM file
			exec " samtools index ${FINAL_BAM_FILE} ${FINAL_BAM_INDEX_FILE}"
			exec "samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}"

			//# =============================
			//# Compute library complexity
			//# =============================
			//# sort by position and strand
			//# Obtain unique count statistics

			PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

			//# PBC File output
			//# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair

			exec """
				bamToBed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf “%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}
			"""
			exec "rm ${FILT_BAM_FILE}"

		}
	}

}

bam_to_tagalign = {

	doc "bam_to_tagalign: ${OFPREFIX}"

	if ( PAIRED_END )
		branch.FINAL_TA_FILE="${FINAL_BAM_PREFIX}.PE2SE.tagAlign.gz"
	else
		branch.FINAL_TA_FILE="${FINAL_BAM_PREFIX}.SE.tagAlign.gz"

	branch.FINAL_BEDPE_FILE="${FINAL_NMSRT_BAM_PREFIX}.bedpe.gz"
	branch.SUBSAMPLED_TA_FILE="${OFPREFIX}.filt.nodup.sample.${NREADS_PER_MILLION}.MATE1.tagAlign.gz"

	if ( PAIRED_END ) {

		produce( FINAL_TA_FILE, FINAL_BEDPE_FILE, SUBSAMPLED_TA_FILE ) {
			//# ===================
			//# Create tagAlign file
			//# ===================
	
			exec """
				bamToBed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{\$4="N";\$5="1000";print \$0}' | gzip -c > ${FINAL_TA_FILE}
			"""
			//# ================
			//# Create BEDPE file
			//# ================
			exec "bamToBed -bedpe -mate1 -i ${FINAL_NMSRT_BAM_FILE} | gzip -c > ${FINAL_BEDPE_FILE}"

			//# =================================
			//# Subsample tagAlign file
			//# Restrict to one read end per pair for CC analysis
			//# ================================
			exec """
				zcat ${FINAL_BEDPE_FILE} | grep -v "chrM" | shuf -n ${NREADS} | awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$3,"N","1000",\$9}' | gzip -c > ${SUBSAMPLED_TA_FILE}
			"""
		}
	}
	else {
		produce( FINAL_TA_FILE, SUBSAMPLED_TA_FILE ) {
			exec """
				bamToBed -i ${FINAL_BAM_FILE} | awk 'BEGIN{OFS="\t"}{\$4="N";\$5="1000";print \$0}' | gzip -c > ${FINAL_TA_FILE}
			"""

			# =================================
			# Subsample tagAlign file
			# ================================
			exec """
				zcat ${FINAL_TA_FILE} | grep -v "chrM" | shuf -n ${NREADS} | gzip -c > ${SUBSAMPLED_TA_FILE}
			"""
		}
	}

}

xcor = {
	doc "xcor: ${OFPREFIX}"

	branch.CC_SCORES_FILE="${SUBSAMPLED_TA_FILE}.cc.qc"
	branch.CC_PLOT_FILE="${SUBSAMPLED_TA_FILE}.cc.plot.pdf"

	produce( CC_SCORES_FILE ) {

		//# CC_SCORE FILE format
		//# Filename <tab> numReads <tab> estFragLen <tab> corr_estFragLen <tab> PhantomPeak <tab> corr_phantomPeak <tab> argmin_corr <tab> min_corr <tab> phantomPeakCoef <tab> relPhantomPeakCoef <tab> QualityTag

//		exec "${RSCRIPT} \$(which run_spp.R) -c=${SUBSAMPLED_TA_FILE} -p=${NTHREADS} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE}"
		exec "${RSCRIPT} \$(which run_spp.R) -c=${SUBSAMPLED_TA_FILE} -filtchr=chrM -savp=${CC_PLOT_FILE} -out=${CC_SCORES_FILE}"

		exec "sed -r 's/,[^\t]+//g' ${CC_SCORES_FILE} > tmp.${OFPREFIX}"
		exec "mv tmp.${OFPREFIX} ${CC_SCORES_FILE}"
	}
}

self_pseudo_rep = {
	doc "self_psedorep: ${OFPREFIX}"

	branch.PR_PREFIX="${OFPREFIX}.filt.nodup"
	branch.PR1_TA_FILE="${PR_PREFIX}.PE2SE.pr1.tagAlign.gz"
	branch.PR2_TA_FILE="${PR_PREFIX}.PE2SE.pr2.tagAlign.gz"

	produce( PR1_TA_FILE, PR2_TA_FILE  ) {
		//# ========================
		//# Create pseudoReplicates
		//# =======================
	
		//# Get total number of read pairs
		//# Shuffle and split BEDPE file into 2 equal parts
		exec """
			nlines=\$( zcat ${FINAL_BEDPE_FILE} | wc -l );
			nlines=\$(( (nlines + 1) / 2 ));
			zcat ${FINAL_BEDPE_FILE} | shuf | split -d -l ${nlines} - ${PR_PREFIX} # Will produce ${PR_PREFIX}00 and ${PR_PREFIX}01"
		"""
	
		//# Convert read pairs to reads into standard tagAlign file
		exec """
			awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' "${PR_PREFIX}00" | gzip -c > ${PR1_TA_FILE}
		"""
		exec """
			rm "${PR_PREFIX}00"
		"""
		exec """
			awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",\$1,\$2,\$3,\$9,\$4,\$5,\$6,\$10}' "${PR_PREFIX}01" | gzip -c > ${PR2_TA_FILE}
		"""
		exec """
			rm "${PR_PREFIX}00"
		"""
	}
}

pool = {
	var IN : "Ctl0"

	from("*${IN}*") {
	
		branch.DATASET_PREFIX = input1.replaceAll(".filt.nodup.PE2SE.pr1.tagAlign.gz","").replaceAll("_Ctl0","").replaceAll("_Ctl1","").replaceAll("_Rep1","").replaceAll("_Rep2","")
	
		doc "pool: ${DATASET_PREFIX}_${IN}"
		//ENCFF002ELJENCFF002ELK_Ctl0_Rep1.filt.srt.nodup.PE2SE.tagAlign.gz
	
		branch.REP1_TA_FILE=  "${DATASET_PREFIX}_${IN}_Rep1.filt.srt.nodup.PE2SE.tagAlign.gz"
		branch.REP2_TA_FILE=  "${DATASET_PREFIX}_${IN}_Rep2.filt.srt.nodup.PE2SE.tagAlign.gz"

		branch.POOLED_TA_FILE="${DATASET_PREFIX}_${IN}_Rep0.filt.srt.nodup.PE2SE.tagAlign.gz"
	
		branch.REP1_PR1_TA_FILE="${DATASET_PREFIX}_${IN}_Rep1.filt.nodup.PE2SE.pr1.tagAlign.gz"
		branch.REP1_PR2_TA_FILE="${DATASET_PREFIX}_${IN}_Rep1.filt.nodup.PE2SE.pr2.tagAlign.gz"
	
		branch.REP2_PR1_TA_FILE="${DATASET_PREFIX}_${IN}_Rep2.filt.nodup.PE2SE.pr1.tagAlign.gz"
		branch.REP2_PR2_TA_FILE="${DATASET_PREFIX}_${IN}_Rep2.filt.nodup.PE2SE.pr2.tagAlign.gz"

		branch.PPR1_TA_FILE="${DATASET_PREFIX}_${IN}_Rep0.filt.nodup.PE2SE.pr1.tagAlign.gz"
		branch.PPR2_TA_FILE="${DATASET_PREFIX}_${IN}_Rep0.filt.nodup.PE2SE.pr2.tagAlign.gz"
	
		produce( POOLED_TA_FILE, PPR1_TA_FILE, PPR2_TA_FILE ) {
			//# ========================
			//# Create pooled datasets
			//# =======================
			exec "zcat ${REP1_TA_FILE} ${REP2_TA_FILE} | gzip -c > ${POOLED_TA_FILE}"
			
			//# ========================
			//# Create pooled pseudoreplicates
			//# =======================
			exec "zcat ${REP1_PR1_TA_FILE} ${REP2_PR1_TA_FILE} | gzip -c > ${PPR1_TA_FILE}"
			exec "zcat ${REP1_PR2_TA_FILE} ${REP2_PR2_TA_FILE} | gzip -c > ${PPR2_TA_FILE}"
		}
	}
}

spp = {
	branch.DATASET_PREFIX = input1.replaceAll(".filt.srt.nodup.PE2SE.tagAlign.gz","").replaceAll("_Ctl0","").replaceAll("_Ctl1","").replaceAll("_Rep1","").replaceAll("_Rep2","").replaceAll("_Rep0","")

	doc "spp: ${DATASET_PREFIX}"

	branch.REP1_TA_FILE=  "${DATASET_PREFIX}_Ctl0_Rep1.filt.srt.nodup.PE2SE.tagAlign.gz"
	branch.REP2_TA_FILE=  "${DATASET_PREFIX}_Ctl0_Rep2.filt.srt.nodup.PE2SE.tagAlign.gz"
	branch.POOLED_TA_FILE="${DATASET_PREFIX}_Ctl0_Rep0.filt.srt.nodup.PE2SE.tagAlign.gz"

	branch.CTL_REP1_TA_FILE=  "${DATASET_PREFIX}_Ctl1_Rep1.filt.srt.nodup.PE2SE.tagAlign.gz"
	branch.CTL_REP2_TA_FILE=  "${DATASET_PREFIX}_Ctl1_Rep2.filt.srt.nodup.PE2SE.tagAlign.gz"
	branch.CTL_POOLED_TA_FILE="${DATASET_PREFIX}_Ctl1_Rep0.filt.srt.nodup.PE2SE.tagAlign.gz"

	branch.REP1_NPEAK="${REP1_TA_FILE}.narrowPeak"
	branch.REP1_CC_SCORES_FILE="${REP1_TA_FILE}.cc"
	branch.REP1_CC_PLOT_FILE="${REP1_TA_FILE}.cc.plot.pdf"

	branch.REP2_NPEAK="${REP2_TA_FILE}.narrowPeak"
	branch.REP2_CC_SCORES_FILE="${REP2_TA_FILE}.cc"
	branch.REP2_CC_PLOT_FILE="${REP2_TA_FILE}.cc.plot.pdf"

	branch.POOLED_NPEAK="${POOLED_TA_FILE}.narrowPeak"
	branch.POOLED_CC_SCORES_FILE="${POOLED_TA_FILE}.cc"
	branch.POOLED_CC_PLOT_FILE="${POOLED_TA_FILE}.cc.plot.pdf"

	produce( REP1_NPEAK, REP1_CC_SCORES_FILE, REP1_CC_PLOT_FILE, 
		 REP2_NPEAK, REP2_CC_SCORES_FILE, REP2_CC_PLOT_FILE, 
		 POOLED_NPEAK, POOLED_SCORES_FILE, POOLED_PLOT_FILE ) {

//		multi	"${RSCRIPT} \$(which run_spp.R) -c=${REP1_TA_FILE} -i=${CTL_REP1_TA_FILE} -filtchr=chrM -p=${NTHREADS} -npeak=${NPEAK} -speak=${SPEAK} -savr=${REP1_NPEAK} -savp=${REP1_CC_PLOT_FILE} -out=${REP1_CC_SCORES_FILE}",
//		 	"${RSCRIPT} \$(which run_spp.R) -c=${REP2_TA_FILE} -i=${CTL_REP2_TA_FILE} -filtchr=chrM -p=${NTHREADS} -npeak=${NPEAK} -speak=${SPEAK} -savr=${REP2_NPEAK} -savp=${REP2_CC_PLOT_FILE} -out=${REP2_CC_SCORES_FILE}",
//			"${RSCRIPT} \$(which run_spp.R) -c=${POOLED_TA_FILE} -i=${CTL_POOLED_TA_FILE} -filtchr=chrM -p=${NTHREADS} -npeak=${NPEAK} -speak=${SPEAK} -savr=${POOLED_NPEAK} -savp=${POOLED_CC_PLOT_FILE} -out=${POOLED_CC_SCORES_FILE}"


		multi	"${RSCRIPT} \$(which run_spp.R) -c=${REP1_TA_FILE} -i=${CTL_REP1_TA_FILE} -filtchr=chrM -npeak=${NPEAK} -speak=${SPEAK} -savr=${REP1_NPEAK} -savp=${REP1_CC_PLOT_FILE} -out=${REP1_CC_SCORES_FILE}",
		 	"${RSCRIPT} \$(which run_spp.R) -c=${REP2_TA_FILE} -i=${CTL_REP2_TA_FILE} -filtchr=chrM -npeak=${NPEAK} -speak=${SPEAK} -savr=${REP2_NPEAK} -savp=${REP2_CC_PLOT_FILE} -out=${REP2_CC_SCORES_FILE}",
			"${RSCRIPT} \$(which run_spp.R) -c=${POOLED_TA_FILE} -i=${CTL_POOLED_TA_FILE} -filtchr=chrM -npeak=${NPEAK} -speak=${SPEAK} -savr=${POOLED_NPEAK} -savp=${POOLED_CC_PLOT_FILE} -out=${POOLED_CC_SCORES_FILE}"


	}
}

idr = {
}

run {
	prepare +  // bpipe bug

	[	"%_Ctl0_Rep1.PE*.fastq.gz" * [ align_bwa + post_align_filt + bam_to_tagalign + xcor + self_pseudo_rep ],
		"%_Ctl0_Rep2.PE*.fastq.gz" * [ align_bwa + post_align_filt + bam_to_tagalign + xcor + self_pseudo_rep ],
		"%_Ctl1_Rep1.PE*.fastq.gz" * [ align_bwa + post_align_filt + bam_to_tagalign + xcor + self_pseudo_rep ],
		"%_Ctl1_Rep2.PE*.fastq.gz" * [ align_bwa + post_align_filt + bam_to_tagalign + xcor + self_pseudo_rep ],
	] + 

	[	pool.using( IN : "Ctl0" ),
		pool.using( IN : "Ctl1" ),
	] 
}




