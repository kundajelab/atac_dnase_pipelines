#!/usr/bin/env python
'''
Created on June 11,2012

@author: ferhat
'''
### import statements ###
import sys
import os
import math
import time
import numpy as np
from scipy import *
from scipy.interpolate import Rbf, UnivariateSpline
from scipy import optimize
from optparse import OptionParser
import scipy.special as scsp
import bisect
# R dependencies
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
fdrtool = importr('fdrtool')
import gzip


from pylab import *
from random import *
from scipy.stats.mstats import mquantiles
import myStats
import myUtils

##### global variables shared by functions ######
chrList=[] # list of all chromosomes
listOfMappableFrags=[] # list of all mappable fragments
# dictkey should be chr-locus1-locus2 and the return values are [interactionDistance, interactionCount]
possiblePairsPerDistance={} # dictionary of possible intra-chr pairs

# these are going to be the important numbers that we compute
possibleInterAllCount=0
observedInterAllCount=0
observedInterAllSum=0
# intra chromosomal interactions in range
possibleIntraInRangeCount=0
observedIntraInRangeCount=0
observedIntraInRangeSum=0
# all possible intra chromosomal interactions
possibleIntraAllCount=0
observedIntraAllCount=0
observedIntraAllSum=0

baselineIntraChrProb=0	# 1.0/possibleIntraAllCount
baselineInterChrProb=0	# 1.0/possibleInterAllCount

minObservedGenomicDist=500000000 # some number bigger than the biggest chromosome length
maxObservedGenomicDist=0 

#distScaling just avoids overflow - but is necessary for large genomes
distScaling=10000.0 
toKb=10**-3
toMb=10**-6
toProb=10**5

### Parameters that can be played with for the SPLINE FIT
overSample=5 # can be changed to have more/less overfitted splines
####
#########################
versionStr="fit-hi-c version 1.0.1. \nA tool for assigning statistical confidence estimates to intra-chromosomal \ncontact maps produced by genome architecture assays. \n\nReleased on January 19, 2014. \nMethod developed by Ferhat Ay, Timothy Bailey and William Noble. \nImplemented by Ferhat Ay (ferhatay@uw.edu). \n\nCopyright (c), 2012, University of Washington. \nThis software is offered under an MIT license. \nFor details: http://opensource.org/licenses/MIT\n"


def main():
	### parse the command line arguments
	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	parser.add_option("-f", "--fragments", dest="fragsfile",
					  help="midpoints (or start indices) of the fragments are read from FRAGSFILE")
	parser.add_option("-i", "--interactions", dest="intersfile",
					  help="interactions between fragment pairs are read from INTERSFILE")
	parser.add_option("-o", "--outdir", dest="outdir",
					  help="where the output files will be written")
	parser.add_option("-t", "--biases", dest="biasfile",
					  help="OPTIONAL: biases calculated by ICE for each locus are read from BIASFILE")
	parser.add_option("-p", "--passes", dest="noOfPasses",type="int",
					  help="OPTIONAL: number of passes after the initial (before) fit. DEFAULT is 1 (after)")
	parser.add_option("-b", "--noOfBins", dest="noOfBins", type="int",
					  help="OPTIONAL: number of equal-occupancy (count) bins. Default is 100")
	parser.add_option("-m", "--mappabilityThres", dest="mappabilityThreshold", type="int",
					  help="OPTIONAL: minimum number of hits per locus that has to exist to call it mappable. DEFAULT is 1.")
	parser.add_option("-l", "--lib", dest="libname",
					  help="OPTIONAL: Name of the library that is analyzed to be used for plots.")
	parser.add_option("-U", "--upperbound", dest="distUpThres", type="int",
					  help="OPTIONAL: upper bound on the intra-chromosomal distance range (unit: base pairs). DEFAULT no limit.")
	parser.add_option("-L", "--lowerbound", dest="distLowThres", type="int",
					  help="OPTIONAL: lower bound on the intra-chromosomal distance range (unit: base pairs). DEFAULT no limit.")
	parser.add_option("-v", "--visual",
					  action="store_true", dest="visual", help="OPTIONAL: use this flag for generating plots. DEFAULT is False.")
	parser.add_option("-q", "--quiet",
					  action="store_false", dest="visual", help="OPTIONAL: use this flag for omitting plots. DEFAULT behavior." )
	parser.add_option("-V", "--version", action="store_true", dest="version", 
					  help=versionStr)
	parser.set_defaults(visual=False, noOfBins=100, distLowThres=-1, distUpThres=-1, mappabilityThreshold=1,noOfPasses=1,
	discBinsize=5000,libname="",biasfile='none', version=False)
	(options, args) = parser.parse_args()
	if len(args) != 0:
		parser.error("incorrect number of arguments")
	if options.version==True:
		print versionStr
		return

	# Set to True for generating plots and False for omitting them
	if options.visual==True:
		# imports related to matplotlib to generate plots
		import matplotlib
		#matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
		#### matplotlib fontsize settings
		plt.rcParams['font.size']=17
		plt.rcParams['axes.labelsize']='x-large'
		plt.rcParams['xtick.labelsize']='large'
		plt.rcParams['ytick.labelsize']='large'
		plt.rcParams['figure.subplot.hspace']=0.4
		plt.rcParams['figure.subplot.bottom']=0.12
		plt.rcParams['figure.subplot.left']=0.15
		plt.rcParams['figure.subplot.right']=0.94
		plt.rcParams['figure.subplot.top']=0.92
	##

	global noOfBins
	global distUpThres
	global distLowThres
	global useBinning
	global libname
	global mappabilityThreshold
	global noOfPasses
	global useInters
	global discBinsize
	global residualFactor
	global outdir
	global visual 
	noOfBins=options.noOfBins # 100 by default 
	distUpThres=options.distUpThres # -1 by default, means no upper bound
	distLowThres=options.distLowThres # -1 by default, means no lower bound
	mappabilityThreshold=options.mappabilityThreshold # 1 by default
	useBinning=True # This is no more an option.
	useInters=False # This is no more an option.
	libname=options.libname
	noOfPasses=options.noOfPasses
	discBinsize=options.discBinsize
	compareMethods=False # This is no more an option.
	residualFactor=-1 # This is no more an option.
	outdir=options.outdir
	visual=options.visual

	# read the mandatory input files -f and -i
	generate_FragPairs(options.fragsfile)
	biasDic={}
	if options.biasfile!='none':
		biasDic=read_ICE_biases(options.biasfile)
	sortedInteractions=read_All_Interactions(options.intersfile,biasDic)

	### DO THE FIRST PASS ###
	# calculate priors using original fit-hic and plot with standard errors
	x,y,yerr=calculate_Probabilities(sortedInteractions,[0 for i in range(len(sortedInteractions))],libname+".fithic_pass1")
	# now fit spline to the data 
	splineXinit,splineYinit,splineResidual,isOutlier,splineFDRxinit,splineFDRyinit=fit_Spline(x,y,yerr,options.intersfile,sortedInteractions,biasDic,libname+".spline_pass1",1)

	### DO THE NEXT PASSES IF REQUESTED ###
	for i in range(2,2+noOfPasses):
		x,y,yerr=calculate_Probabilities(sortedInteractions,isOutlier,libname+".fithic_pass"+repr(i))
		splineX,splineY,splineResidual,isOutlier,splineFDRx,splineFDRy=fit_Spline(x,y,yerr,options.intersfile,sortedInteractions,biasDic,libname+".spline_pass"+repr(i),i)

	print
	sys.stderr.write("\nExecution of fit-hic completed successfully. \n\n") 
	return # from main

def read_ICE_biases(infilename):
	sys.stderr.write("\n\nReading ICE biases. \n")
	biasDic={}
	
	rawBiases=[]
	infile =gzip.open(infilename, 'r')
	for line in infile:
		words=line.rstrip().split()
		chr=words[0]; midPoint=int(words[1]); bias=float(words[2])
		if bias!=1.0:
			rawBiases.append(bias)
	infile.close()
	#sys.stderr.write("\n\nReading ICE biases. \n")
	botQ,med,topQ=mquantiles(rawBiases,prob=[0.05,0.5,0.95])
	#sys.stderr.write("5th quantile of biases: "+str(botQ)+"\n")
	#sys.stderr.write("50th quantile of biases: "+str(med)+"\n")
	#sys.stderr.write("95th quantile of biases: "+str(topQ)+"\n")

	#normFactor=sum(rawBiases)/len(rawBiases)
	infile =gzip.open(infilename, 'r')
	totalC=0
	discardC=0
	for line in infile:
		words=line.rstrip().split()
		chr=words[0]; midPoint=int(words[1]); bias=float(words[2])
		# extra conditions
		#if bias<(botQ/2.0):
		if bias<0.5:
			bias=-1 #botQ
			discardC+=1
		elif bias>2:
			bias=-1 #topQ
			discardC+=1
		#
		totalC+=1
		if chr not in biasDic:
			biasDic[chr]={}
		if midPoint not in biasDic[chr]:
			biasDic[chr][midPoint]=bias
	infile.close()
	sys.stderr.write("Out of " + str(totalC) + " loci " +str(discardC) +" were discarded with biases not in range [0.5 2]\n\n" )

	return biasDic # from read_ICE_biases

def fit_Spline(x,y,yerr,infilename,sortedInteractions,biasDic,figname,passNo):
	sys.stderr.write("\nFit a univariate spline to the probability means\n")
	sys.stderr.write("------------------------------------------------------------------------------------\n")
	sys.stderr.write("baseline intra-chr probability: " + repr(baselineIntraChrProb)+ "\tbaseline inter-chr probability: " + repr(baselineInterChrProb)+"\n")
	# xi and yi will be used only for visualization purposes
	# acutal fit and residual is all done on vectors x and y
	xi = np.linspace(min(x), max(x), overSample*len(x))

	# assume residualFactor==-1: 
	splineError=min(y)*min(y)

	# use fitpack2 method -fit on the real x and y from equal occupancy binning
	ius = UnivariateSpline(x, y, s=splineError)
	yi = ius(xi)

	#### POST-PROCESS THE SPLINE TO MAKE SURE IT'S NON-INCREASING
	### NOW I DO THIS BY CALLING AN R function CALLED MONOREG 
	### This does the isotonic regression using option antitonic to make sure 
	### I get monotonically decreasing probabilites with increasion genomic distance 

	tempMaxX=max(x)
	tempMinX=min(x)
	tempList=sorted(list(set([int(i[0]) for i in sortedInteractions])))
	splineX=[]
	### The below for loop will make sure nothing is out of range of [min(x) max(x)]
	### Therefore everything will be within the range where the spline is defined
	for i in tempList:
		if tempMinX<=i and i<=tempMaxX:
			splineX.append(i)
	# END for
	#print len(splineX)
	splineY=ius(splineX)

	# R vector format
	rSplineX=ro.FloatVector(splineX)
	rSplineY=ro.FloatVector(splineY)
	rMonoReg=ro.r['monoreg']
	# do the antitonic regression
	allRres=rMonoReg(rSplineX,rSplineY,type="antitonic")
	rNewSplineY=allRres[3]
	# convert data back to Python format
	newSplineY=[]
	diff=[]
	diffX=[]
	for i in range(len(rNewSplineY)):
		newSplineY.append(rNewSplineY[i])
		if (splineY[i]-newSplineY[i]) > 0:
			diff.append(splineY[i]-newSplineY[i])
			diffX.append(splineX[i])
	# END for
	#print len(splineX)
	
	residual =sum([i*i for i in (y - ius(x))])

	if visual==True:
		### Now plot the results
		sys.stderr.write("Plotting %s" % figname + ".png\n")
		plt.clf()
		fig = plt.figure()
		ax = fig.add_subplot(2,1,1)
		plt.plot(myUtils.scale_a_list(splineX,toKb), myUtils.scale_a_list(newSplineY,toProb),'g-',label="spline-"+str(passNo),linewidth=2)
		plt.errorbar(myUtils.scale_a_list(x,toKb),myUtils.scale_a_list(y,toProb),myUtils.scale_a_list(yerr,toProb),fmt='r.',label="Mean with std. error",linewidth=2) 

		if useInters:
			plt.plot(myUtils.scale_a_list(x,toKb),myUtils.scale_a_list([baselineIntraChrProb for i in x],toProb),'k-',label="Baseline intra-chromosomal")
			plt.plot(myUtils.scale_a_list(x,toKb),myUtils.scale_a_list([baselineIntraChrProb for i in x],toProb),'b-',label="Baseline inter-chromosomal")
		plt.ylabel('Contact probability (x10$^{-5}$)',fontsize='large')
		plt.xlabel('Genomic distance (kb)',fontsize='large')
		if distLowThres>-1 and distUpThres>-1:
			plt.xlim(myUtils.scale_a_list([distLowThres, distUpThres],toKb))
		plt.gca().yaxis.set_major_locator( MaxNLocator(nbins = 3, prune=None))
		ax.legend(loc="upper right")

		ax = fig.add_subplot(2,1,2)

		plt.loglog(splineX,newSplineY,'g-')
		plt.errorbar(x, y, yerr=yerr, fmt='r.') # Data
		if useInters:
			plt.loglog(x,[baselineIntraChrProb for i in x],'k-')
			plt.loglog(x,[baselineIntraChrProb for i in x],'b-')
		if distLowThres>-1 and distUpThres>-1:
			plt.xlim([distLowThres, distUpThres])
		plt.ylabel('Contact probability (log-scale)',fontsize='large')
		plt.xlabel('Genomic distance (log-scale)',fontsize='large')

		plt.savefig(outdir+'/'+figname+'.png')

	# NOW write the calculated pvalues and corrected pvalues in a file 
	infile =gzip.open(infilename, 'r')
	intraInRangeCount=0
	intraOutOfRangeCount=0
	intraVeryProximalCount=0
	interCount=0
	sys.stderr.write("distLowThres " + repr(distLowThres) + "\tdistUpThres " + repr(distUpThres) +"\n")
	p_vals=[]
	q_vals=[]
	for line in infile:
		words=line.rstrip().split()
		interxn=myUtils.Interaction([words[0], int(words[1]), words[2], int(words[3])])
		interxn.setCount(int(words[4]))
		chr1=words[0]
		chr2=words[2]
		midPoint1=int(words[1])
		midPoint2=int(words[3])

		bias1=1.0; bias2=1.0;  # assumes there is no bias to begin with
		# if the biasDic is not null sets the real bias values
		if len(biasDic)>0:
			if chr1 in biasDic and midPoint1 in biasDic[chr1]:
				bias1=biasDic[chr1][midPoint1]
			if chr2 in biasDic and midPoint2 in biasDic[chr2]:
				bias2=biasDic[chr2][midPoint2]

		if (bias1<0 or bias2<0) and interxn.type!='inter':
			prior_p=1.0
			p_val=1.0
			p_vals.append(p_val)
		elif interxn.getType(distLowThres,distUpThres)=='intraInRange': 
			# make sure the interaction distance is covered by the probability bins
			distToLookUp=max(interxn.distance,min(x))
			distToLookUp=min(distToLookUp,max(x))
			i=min(bisect.bisect_left(splineX, distToLookUp),len(splineX)-1) 
			#prior_p=newSplineY[i]
			prior_p=newSplineY[i]*(bias1*bias2) # biases added in the picture
			intraInRangeCount +=1
			############# THIS HAS TO BE interactionCount-1 ##################
			p_val=scsp.bdtrc(interxn.hitCount-1,observedIntraInRangeSum,prior_p)
			p_vals.append(p_val)

		elif interxn.getType(distLowThres,distUpThres)=='intraShort':
			prior_p=1.0
			p_val=1.0
			intraVeryProximalCount +=1
			p_vals.append(p_val)

		elif interxn.getType(distLowThres,distUpThres)=='intraLong':
			# out of range bigger than distUpThres
			# use the prior of the baseline intra-chr interaction probability
			prior_p=1.0 #baselineIntraChrProb*(bias1*bias2)  # biases added in the picture
			p_val=scsp.bdtrc(interxn.hitCount-1,observedIntraAllSum,prior_p)
			intraOutOfRangeCount +=1
			p_vals.append(p_val)

		else:
			if useInters:
				#prior_p=baselineIntraChrProb
				prior_p=baselineInterChrProb*(bias1*bias2) # biases added in the picture
				############# THIS HAS TO BE interactionCount-1 ##################
				p_val=scsp.bdtrc(interxn.hitCount-1,observedInterAllSum,prior_p)
				interCount +=1
				p_vals.append(p_val)
	# END for
	infile.close()

	# Do the BH FDR correction 
	if useInters:
		q_vals=myStats.benjamini_hochberg_correction(p_vals, possibleInterAllCount+possibleIntraAllCount)
		sys.stderr.write("possibleInterAllCount+possibleIntraAllCount " + repr(possibleInterAllCount+possibleIntraAllCount)+"\n")
	else:
		q_vals=myStats.benjamini_hochberg_correction(p_vals, possibleIntraInRangeCount)
		sys.stderr.write("possibleIntraInRangeCount " + repr(possibleIntraInRangeCount)+"\n")

	infile =gzip.open(infilename, 'r')
	outfile =gzip.open(outdir+'/'+figname+'.significances.txt.gz', 'w')
	sys.stderr.write("Writing p-values to file %s" % figname + ".significances.txt.gz\n")
	count=0
	outfile.write("chr1\tfragmentMid1\tchr2\tfragmentMid2\tcontactCount\tp-value\tq-value\n")

	for line in infile:
		words=line.rstrip().split()
		chrNo1=words[0]
		midPoint1=int(words[1])
		chrNo2=words[2]
		midPoint2=int(words[3])
		interactionCount=int(words[4])
		p_val=p_vals[count]
		q_val=q_vals[count]
		
		if useInters==False and chrNo1==chrNo2: # intra
			interactionDistance=abs(midPoint1-midPoint2) # dist 
			if myUtils.in_range_check(interactionDistance,distLowThres,distUpThres):
				outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\n" % (str(chrNo1),midPoint1,str(chrNo2),midPoint2,interactionCount,p_val,q_val))
		elif useInters==True and chrNo1!=chrNo2:
			outfile.write("%s\t%d\t%s\t%d\t%d\t%e\t%e\n" % (str(chrNo1),midPoint1,str(chrNo2),midPoint2,interactionCount,p_val,q_val))
		#outfile.write("ALL\t%s\t%d\t%s\t%d\t%d\t%e\t%e\n" % (str(chrNo1),midPoint1,str(chrNo2),midPoint2,interactionCount,p_val,q_val))

		count+=1
	# END for - printing pvals and qvals for all the interactions
	outfile.close()

	isOutlier=[]
	distsBelow=[]
	distsAbove=[]
	intcountsBelow=[]
	intcountsAbove=[]
	belowThresCount=0
	aboveThresCount=0
	outlierThres=1.0/possibleIntraInRangeCount
	for interactionDistance,interactionCount,bias12 in sortedInteractions:
		# make sure the interaction distance is covered by the probability bins
		distToLookUp=max(interactionDistance,min(x))
		distToLookUp=min(distToLookUp,max(x))
		i=min(bisect.bisect_left(splineX, distToLookUp),len(splineX)-1) 
		prior_p=newSplineY[i]*float(bias12) # biases added in the picture
		############# THIS HAS TO BE interactionCount-1 ##################
		p_val=scsp.bdtrc(interactionCount-1,observedIntraInRangeSum,prior_p)
		if p_val < outlierThres:
			distsBelow.append(interactionDistance)
			intcountsBelow.append(interactionCount)
			isOutlier.append(1)
			belowThresCount +=1
		else:
			distsAbove.append(interactionDistance)
			intcountsAbove.append(interactionCount)
			isOutlier.append(0)
			aboveThresCount +=1
	# END for - doing the outlier check for all interactions in sortedInteractions


	if visual==True:
		sys.stderr.write("Plotting results of extracting outliers to file %s" % figname + ".extractOutliers.png\n")
		plt.clf()
		fig = plt.figure()
		ax = fig.add_subplot(111)
		downsample=30 # for the non-outliers
		randIndcsAbove=sample([i for i in range(len(intcountsAbove))],len(intcountsAbove)/downsample)
		randIndcsAbove=sorted(randIndcsAbove)
		downsample=20 # for the outliers
		randIndcsBelow=sample([i for i in range(len(intcountsBelow))],len(intcountsBelow)/downsample)
		randIndcsBelow=sorted(randIndcsBelow)

		plt.plot(myUtils.scale_a_list([distsBelow[i] for i in randIndcsBelow],toKb),[intcountsBelow[i] for i in randIndcsBelow], 'r.',label="Outliers (p-value < 1/M)")
		plt.plot(myUtils.scale_a_list(splineX+[maxObservedGenomicDist],toKb),[newSplineY[i]*observedIntraInRangeSum	for i in range(len(newSplineY))]+[newSplineY[-1]*observedIntraInRangeSum], 'g-', label="spline-"+str(passNo)+" (x N)", linewidth=2.5)

		plt.xlabel('Genomic distance (kb)')
		plt.ylabel('Contact counts')
		print(repr(len(intcountsBelow))+"\t"),
		## this limits y-axis of the hit count plots
		if len(intcountsBelow)>0:
			plt.ylim([0,min(max(intcountsBelow),1500)])
		if distLowThres>-1 and distUpThres>-1:
			plt.xlim([0, distUpThres*toKb])
		ax.legend(loc="upper right",fancybox=True)
		plt.savefig(outdir+'/'+figname+'.extractOutliers.png')

	sys.stderr.write("intraInRangeCount " + repr(intraInRangeCount)+"\tintraOutOfRangeCount " +\
		repr(intraOutOfRangeCount)+"\tintraVeryProximalCount " + repr(intraVeryProximalCount) +"\tinterCount " + repr(interCount)+"\n")

	if visual==True:
		sys.stderr.write("Plotting q-values to file %s" % figname + ".qplot.png\n")
	minFDR=0.0
	maxFDR=0.05
	increment=0.001
	FDRx,FDRy=plot_qvalues(q_vals,minFDR,maxFDR,increment,figname+".qplot")

	infile.close()

	return [splineX, newSplineY, residual, isOutlier, FDRx, FDRy] # from fit_Spline


def calculate_Probabilities(sortedInteractions,isOutlier,figname):

	sys.stderr.write("\nCalculating probability means and standard deviations by equal occupancy binning of interaction data\n")
	sys.stderr.write("------------------------------------------------------------------------------------\n")
	
	outfile =open(outdir+'/'+figname+'.txt', 'w')
	
	# total interaction count to put on top of the plot
	# this may be different than observedIntraInRangeSum for the second iteration of fit-hic
	totalInteractionCountForPlot=0
	lcount=0
	for eachrow in sortedInteractions:
		if isOutlier[lcount]==0:
			totalInteractionCountForPlot += eachrow[1]
		lcount+=1
	# END for
	desiredPerBin=(observedIntraInRangeSum)/noOfBins
	sys.stderr.write("observedIntraInRangeSum\t"+repr(observedIntraInRangeSum)+ "\tdesiredPerBin\t" +repr(desiredPerBin)+"\tnoOfBins\t"+repr(noOfBins)+"\n")

	# the following five lists will be the print outputs
	x=[] # avg genomic distances of bins
	y=[] # avg interaction probabilities of bins
	yerr=[] # stderrs of bins
	pairCounts=[] # number of pairs in bins
	interactionTotals=[] # number of interactions (reads) in bins

	# the following variables will be used to calculate the above five lists
	noOfPairsForBin=0
	meanCountPerPair=0
	M2=0
	interactionTotalForBin=0
	interactionTotalForBinTermination=0
	distanceTotalForBin=0
	lastDistanceForBin=-1
	lastInteraction=lcount 
	lcount=0 # this will increase by eachrow in sortedInteractions

	for eachrow in sortedInteractions:
		interactionDistance=eachrow[0]
		interactionCount=eachrow[1]
 
 		# if one bin is full or it's the last bin
		if noOfPairsForBin>0 and ((useBinning==False and lastDistanceForBin!=-1 and lastDistanceForBin!=interactionDistance) or\
			(useBinning==True and lastDistanceForBin!=-1 and interactionTotalForBinTermination >= desiredPerBin and\
			lastDistanceForBin!=interactionDistance) or lcount==lastInteraction): 

			# calculate the things that need to be calculated
			avgDistance=(distanceTotalForBin/noOfPairsForBin)*distScaling
			meanProbabilityObsv=(meanCountPerPair*1.0)/observedIntraInRangeSum
			se_p=meanProbabilityObsv
			# update se_p if there are more than 1 pairs in the bin
			if noOfPairsForBin>1:
				var=M2/(noOfPairsForBin-1)
				sd=math.sqrt(var)
				se=sd/math.sqrt(noOfPairsForBin)
				se_p=se/observedIntraInRangeSum
			# END if

			# append the calculated vals to corresponding lists
			x.append(float(avgDistance))
			y.append(float(meanProbabilityObsv))
			yerr.append(float(se_p))
			pairCounts.append(noOfPairsForBin)
			interactionTotals.append(interactionTotalForBin)
	
			# now that we saved what we need
			# set the values back to defaults and go on to the next bin
			noOfPairsForBin=0
			meanCountPerPair=0
			M2=0
			interactionTotalForBin=0
			interactionTotalForBinTermination=0
			distanceTotalForBin=0
			lastDistanceForBin=-1
		# END if - that checks whether the bin is full etc.

		# Now go back to processing the read values of interactionDistance and interactionCount
		# this check is necessary for the second pass of fit-hic
		# we want to only use the non-outlier interactions in our
		# probability calculation
		if isOutlier[lcount]==0:
			distanceTotalForBin +=interactionDistance/distScaling
			interactionTotalForBin +=interactionCount
			noOfPairsForBin +=1
			delta=interactionCount-meanCountPerPair
			meanCountPerPair += (delta*1.0) / noOfPairsForBin
			M2 +=delta*(interactionCount-meanCountPerPair)
		# END if
		interactionTotalForBinTermination +=interactionCount
		lcount +=1
		lastDistanceForBin=interactionDistance
	# END for over sortedInteractions

	if visual==True:
		sys.stderr.write("Plotting %s" % figname + ".png\n")
		plt.clf()
		fig = plt.figure()
		ax = fig.add_subplot(111)
		plt.plot(myUtils.scale_a_list(x,toKb),myUtils.scale_a_list(y,toProb),'ro',label="Mean")
		plt.errorbar(myUtils.scale_a_list(x,toKb),myUtils.scale_a_list(y,toProb),myUtils.scale_a_list(yerr,toProb),fmt='k.', label="Standard error")
		#plt.ylabel('Probability (1e-5)')
		plt.ylabel('Contact probability (x10$^{-5}$)')
		plt.xlabel('Genomic distance (kb)')
		titleStr='Binning observed interactions using equal occupancy bins.\n No. of bins: '\
			+str(noOfBins) +', Library: ' + str(libname)+ ', No. of interactions: ' +str(observedIntraInRangeSum)
		plt.title(titleStr,size='small')
		ax.legend(loc="upper right")
		plt.savefig(outdir+'/'+figname+'.png')

	sys.stderr.write("Writing %s" % figname + ".txt\n")
	
	outfile.write("avgGenomicDist\tcontactProbability\tstandardError\tnoOfLocusPairs\ttotalOfContactCounts\n")
	for i in range(len(x)):
		outfile.write("%d" % x[i] + "\t"+"%.2e" % y[i]+ "\t" + "%.2e" % yerr[i] + "\t" +"%d" % pairCounts[i] + "\t" +"%d" % interactionTotals[i]+"\n")
	outfile.close()
	return [x,y,yerr] # from calculate_Probabilities

def read_All_Interactions(infilename,biasDic):
	sys.stderr.write("\nReading all the interactions and then sorting the intra chr ones in range according to genomic distance\n")
	sys.stderr.write("------------------------------------------------------------------------------------\n")

	# global variables initialized by this function
	global observedIntraAllSum
	global observedIntraAllCount
	global observedIntraInRangeSum
	global observedIntraInRangeCount
	global observedInterAllSum
	global observedInterAllCount
	global minObservedGenomicDist
	global maxObservedGenomicDist

	#read the interactions file - create a two dimensional numpy array with each row is a [distance,count] pair
	infile =gzip.open(infilename, 'r')
	for line in infile:
		words=line.rstrip().split()
		interxn=myUtils.Interaction([words[0], int(words[1]), words[2], int(words[3])])
		interxn.setCount(int(words[4]))
		chrIndex1=chrList.index(interxn.chr1)
		chrIndex2=chrList.index(interxn.chr2)
		chr1=words[0]
		chr2=words[2]
		midPoint1=int(words[1])
		midPoint2=int(words[3])

		bias1=1.0; bias2=1.0;  # assumes there is no bias to begin with
		# if the biasDic is not null sets the real bias values
		if len(biasDic)>0:
			if chr1 in biasDic and midPoint1 in biasDic[chr1]:
				bias1=biasDic[chr1][midPoint1]
			if chr2 in biasDic and midPoint2 in biasDic[chr2]:
				bias2=biasDic[chr2][midPoint2]


		if interxn.mid1 not in listOfMappableFrags[chrIndex1] or interxn.mid2 not in listOfMappableFrags[chrIndex2]:
			sys.stderr.write("Not-mappable fragment pair: %s\t" %str(interxn.chr1)+"%d\t" % interxn.mid1+ "%s\t" %str(interxn.chr2) +"%d\n" % interxn.mid2)
			continue

		if interxn.type=='inter':
			observedInterAllSum +=interxn.hitCount
			observedInterAllCount +=1
		else: # any type of intra
			observedIntraAllSum +=interxn.hitCount
			observedIntraAllCount +=1
			if interxn.getType(distLowThres,distUpThres)=='intraInRange':
				minObservedGenomicDist=min(minObservedGenomicDist,interxn.distance)
				maxObservedGenomicDist=max(maxObservedGenomicDist,interxn.distance)
				# every pair should already be in the dictionary with a zero interaction count
				dictkey=str(interxn.chr1)+'-'+str(min(interxn.mid1,interxn.mid2))+'-'+str(max(interxn.mid1,interxn.mid2))
				if not dictkey in possiblePairsPerDistance:
					sys.exit("Illegal fragment pair")
				else:
					possiblePairsPerDistance[dictkey]=[interxn.distance,interxn.hitCount,bias1*bias2] #--now with biases
				observedIntraInRangeSum +=interxn.hitCount
				observedIntraInRangeCount +=1
		# END else

	# END for
	infile.close()
	sys.stderr.write("Total of \t"+str(observedIntraAllCount) +" observed intra-chr fragment pairs,\t"\
		+str(observedIntraInRangeCount) +" observed intra-chr fragment pairs in range,\t"\
		+str(observedInterAllCount) +" observed inter-chr fragment pairs\n" )
	sys.stderr.write("Total of \t"+str(observedIntraAllSum) +" observed intra-chr read counts,\t"\
		+str(observedIntraInRangeSum) +" observed intra-chr read counts in range,\t"\
		+str(observedInterAllSum) +" observed inter-chr read counts\n" )
	sys.stderr.write("Range of observed genomic distances	[%d	%d]" % (minObservedGenomicDist,maxObservedGenomicDist) + "\n")

	# sort the interactions if not already sorted
	sortedInteractions=[]
	for i in possiblePairsPerDistance:
		sortedInteractions.append(possiblePairsPerDistance.get(i))
	
	t=time.time()
	myUtils.sort_by_column(sortedInteractions,0) #in-place sorting according to column index 0 (first column)
	sys.stderr.write("Total time for sorting interactions according to genomic distance: %.3f\n" % (time.time()-t))

	return sortedInteractions #from read_All_Interactions

def generate_FragPairs(infilename):
	sys.stderr.write("\nGenerating all possible intra-chromosomal fragment pairs and counting the number of all possible inter-chr fragment pairs\n")
	sys.stderr.write("------------------------------------------------------------------------------------\n")
	global listOfMappableFrags # two dimensional list with all mappable fragment midpoints for each chr
	global chrList # list of all chromosomes (chrno (type=int))
	global possiblePairsPerDistance # all possible intra-chr fragment pairs
	global possibleInterAllCount # count of all possible inter-chr fragment pairs
	global possibleIntraAllCount # count of all possible intra-chr fragment pairs
	global possibleIntraInRangeCount # count of all possible intra-chr fragment pairs in the range we're interested
	global baselineInterChrProb # 1 divided by all possible inter-chr fragment pairs 
	global baselineIntraChrProb #  1 divided by all possible intra-chr fragment pairs

	listOfMappableFrags=[]
	chrList=[]

	#get the name of the first chr
	infile =gzip.open(infilename, 'r')
	line=infile.readline()
	words=line.rstrip().split()
	currChrNo=words[0] #get the name of first chr
	infile.close()

	# read the fragments file 
	fragsPerChr=[] # temporary list that will be added to listOfMappableFrags for each chr
	totalNoOfFrags=0 # total number of all mappable fragments
	infile =gzip.open(infilename, 'r')
	for line in infile:
		words=line.rstrip().split()
		chrNo=words[0] # can be an integer or a string
		#words[1] ignored
		midPoint=int(words[2])
		hitCount=int(words[3])
		# whenever the name of the chromosome changes 
		if currChrNo!=chrNo:
			listOfMappableFrags.append(fragsPerChr)
			totalNoOfFrags += len(fragsPerChr)
			chrList.append(currChrNo)
			currChrNo = chrNo
			fragsPerChr=[]
		# add the mappable midPoints to the temp fragsPerChr
		if hitCount >= mappabilityThreshold:
			fragsPerChr.append(midPoint)
	#END for

	# handle the last chromosome
	listOfMappableFrags.append(fragsPerChr)
	totalNoOfFrags += len(fragsPerChr)
	chrList.append(currChrNo)
	infile.close()
	
	# create all possible frag pairs 
	possibleInterAllCount=0
	possibleIntraInRangeCount=0
	possibleIntraAllCount=0
	for i in chrList:
		countIntraPairs=0
		chrIndex=chrList.index(i) # get the index of chromosome from the chrList 
		fragsPerChr=(listOfMappableFrags[chrIndex])[:] # get the mappable midpoints for that chr
		tempLen=len(fragsPerChr)
		possibleInterAllCount+= (totalNoOfFrags-tempLen)*tempLen
		# iterate over all possible intra-chr pairs to see which ones qualify as a 'possible' pair
		for x in range(tempLen):
			for y in range(x+1,tempLen):
				interactionDistance=abs(fragsPerChr[x]-fragsPerChr[y])
				if myUtils.in_range_check(interactionDistance,distLowThres,distUpThres):
					countIntraPairs +=1
					dictkey=str(i)+'-'+str(min(fragsPerChr[x],fragsPerChr[y]))+'-'+str(max(fragsPerChr[x],fragsPerChr[y]))
					possiblePairsPerDistance[dictkey]=[interactionDistance,0,1.0] # set count to zero for now and bias to 1.0
				possibleIntraAllCount+=1
			#END for
		#END for
		possibleIntraInRangeCount+=countIntraPairs
		sys.stderr.write("Chromosome " +repr(i) +",\t"+str(tempLen) +" mappable fragments, \t"+str(countIntraPairs)\
		+" possible intra-chr fragment pairs in range,\t" + str((totalNoOfFrags-tempLen)*tempLen) +" possible inter-chr fragment pairs\n")
	#END for

	# divide the possibleInterAllCount by 2 so that every inter-chr interaction is counted only once
	possibleInterAllCount=possibleInterAllCount/2
	sys.stderr.write("Total of \t"+str(possibleIntraInRangeCount) +" possible intra-chr fragment pairs in range,\t"\
	+str(possibleIntraAllCount) +" possible intra-chr fragment pairs,\t"\
	+str(possibleInterAllCount) +" possible inter-chr fragment pairs\n")
	# calculate inter-chr probabilities
	if possibleInterAllCount >0:
		baselineInterChrProb=1.0/possibleInterAllCount
	baselineIntraChrProb=1.0/possibleIntraAllCount

	return # from generate_FragPairs



def plot_qvalues(q_values,minFDR,maxFDR,increment,figname):
	qvalTicks=np.arange(minFDR,maxFDR+increment,increment)
	significantTicks=[0 for i in range(len(qvalTicks))]
	qvalBins=[-1 for i in range(len(q_values))]
	for i, q in enumerate(q_values):
		qvalBins[i]=int(math.floor(q/increment))
	
	for i in range(len(qvalBins)):
		if qvalBins[i]>=len(qvalTicks):
			continue
		significantTicks[qvalBins[i]]+=1
	
	# make it cumulative 
	for i in range(1,len(significantTicks)):
		significantTicks[i]=significantTicks[i]+significantTicks[i-1]
	# shift them by 1
	for i in range(1,len(significantTicks)):
		significantTicks[-1*i]=significantTicks[-1*i-1]
	significantTicks[0]=0

	if visual==True:
		plt.clf()
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		plt.plot(qvalTicks,significantTicks, 'b*-')
		plt.xlabel('FDR threshold')
		plt.ylabel('Significant contacts')
		plt.savefig(outdir+'/'+figname+'.png')

	return [qvalTicks,significantTicks]


if __name__ == "__main__":
	main()

