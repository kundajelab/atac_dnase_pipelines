##############################################################################
### To use the functions in this lib simply import this python module using
### import myUtils
### Then you'll able able to call functions with the proper arguments using
### returnVal=myUtils.func1(arg1,arg2)
##############################################################################
##############################################################################
import numpy as np
import math

################### FUNC scale_a_list  #######################################
#### Given a list this functions multiplies it with the given scaling parameter.
##############################################################################
def scale_a_list(somelist, s):
	return [1.0*i*s for i in somelist]

################### FUNC chr_name_conversion  #####################################
#### Given an identifier for the chromosome name (str) or a chromosome number (int) 
###  and an organism this function converts the identifier to the other representation. 
### Example: 

################### FUNC sort_by_column  #####################################
#### Given a list with 1 or more columns this functions sorts it according to 
### the desired column n [0 len(list)). Does this in-place.
##############################################################################
def sort_by_column(somelist, n):
	somelist[:] = [(x[n], x) for x in somelist]
	somelist.sort()
	somelist[:] = [val for (key, val) in somelist]
	return

################### FUNC chr_name_conversion  #####################################
#### Given an identifier for the chromosome name (str) or a chromosome number (int) 
###  and an organism this function converts the identifier to the other representation. 
### Example: 
### converts 'chrX' or 'X' to 23 for human
### converts 23 to 'chrX' 1 to 'chr1' for human
##############################################################################
def chr_name_conversion(chrIn,org):
	if isinstance(chrIn, int): # int to str
		if org=='human':
			if	chrIn<23 and chrIn>0:
				chrOut='chr'+str(chrIn)
			elif chrIn==23:
				chrOut='chrX'
			elif chrIn==24:
				chrOut='chrY'
			else:
				return 'problem'
		elif org=='mouse':
			if	chrIn<20 and chrIn>0:
				chrOut='chr'+str(chrIn)
			elif chrIn==20:
				chrOut='chrX'
			elif chrIn==21:
				chrOut='chrY'
			else:
				return 'problem'
		else:
			chrOut='chr'+str(chrIn)
	else: # str to int
		if 'chr' in chrIn:
			chrIn=chrIn[:3] # cut the 'chr'
		if org=='human':
			if	chrIn=='X':
				chrOut=23
			elif chrIn=='Y':
				chrOut=24
			else:
				chrOut=int(chrIn)
		elif org=='mouse':
			if	chrIn=='X':
				chrOut=20
			elif chrIn=='Y':
				chrOut=21
			else:
				chrOut=int(chrIn)

	return chrOut

################### FUNC in_range_check  #####################################
####  Check whether the given interactionDistance is within the range we are 
###   interested. Should only be used for intra-chromosomal interactions.
##############################################################################
def in_range_check(interactionDistance,distLowThres,distUpThres):
	if (distLowThres==-1 or (distLowThres>-1 and interactionDistance >distLowThres)) and\
		(distUpThres==-1 or (distUpThres>-1 and interactionDistance <= distUpThres)):
		return True
	return False

################### CLASS Interaction ########################################
####  This class is a container for interactions between two loci as observed 
###   in Hi-C datasets. 
##############################################################################
class Interaction:
	chr1='chr'
	mid1=-1
	chr2='chr'
	mid2=-1
	hitCount=0
	distance=-1
	type='null'
	pval=-1.0
	qval=-1.0
	dictkey='null'
	def __init__(self):
		self.data = []

	def __init__(self, locusPair):
		self.chr1=locusPair[0]
		self.mid1=int(locusPair[1])
		self.chr2=locusPair[2]
		self.mid2=int(locusPair[3])
		if self.chr1==self.chr2:
			self.type='intra'
			self.distance=abs(self.mid1-self.mid2)
		else:
			self.type='inter'

	def setCount(self,x):
		self.hitCount=float(x)
	def setType(self,x):
		self.type=str(x)
	def setPval(self,x):
		self.pval=float(x)
	def setQval(self,x):
		self.qval=float(x)
	def getType(self,distLowThres,distUpThres):
		if (distLowThres==-1 or (distLowThres>-1 and self.distance >distLowThres)) and\
			(distUpThres==-1 or (distUpThres>-1 and self.distance <= distUpThres)):
			self.type='intraInRange'
		elif (distLowThres>-1 and self.distance <= distLowThres):
			self.type='intraShort'
		elif (distUpThres>-1 and self.distance > distUpThres):
			self.type='intraLong'
		return self.type

		#if self.chr1==min(self.chr1,self.chr2):
				
################### CLASS Locus ##############################################
####  This class is a container for loci in general and for Hi-C data.
##############################################################################
class Locus:
	chr='chr'
	mid=-1
	#start=-1
	#end=-1
	type='null'
	hitCount=-1
	#avgGCcontent=-1.0
	#avgMappability=-1.0
	def __init__(self):
		self.data = []

	def __init__(self, locus):
		self.chr=locus[0]
		self.mid=int(locus[1])
	def setCount(self,x):
		self.hitCount=float(x)
	def setType(self,x):
		self.type=str(x)


################### FUNC convert_UCSC_to_bed_format  #####################################
#### Given a locus in UCSC format this function converts it to bed format with 3 fields
### chr1:121-21111 --> ['chr1', 121, 21111]
##############################################################################
def convert_UCSC_to_bed_format(l):
	chr=l[:l.find(':')]
	st=int(l[l.find(':')+1:l.find('-')])
	en=int(l[l.find('-')+1:])
	return (chr,st,en)


################### FUNC read_bed_for_chrlist  #####################################
#### Reads a bed file and extracts a dictionary for each chromosome given in the chromosomeList
### <noOfAdditionalFields> dditional fields are going to be considered on top of first 3 fields.
##############################################################################
def read_bed_for_chrlist(bedFile,chromosomeList,noOfAdditionalFields):
	dic={}
	for ch in chromosomeList:
		dic[ch]=[]

	infile=open(bedFile,'r')
	for line in infile:
		## line is too short to be real
		if len(line)<3:
			break
		words=line.rstrip().split()
		thisChr=words[0]
		if thisChr not in chromosomeList:
			continue
		thisStart=int(words[1])
		thisEnd=int(words[2])
		additions=[]
		if noOfAdditionalFields>0:
			for n in range(noOfAdditionalFields):
				additions.append(words[3+n])
		dic[thisChr].append([thisStart, thisEnd]+additions)
	infile.close()
	for ch in chromosomeList:
		dic[ch]=sorted(dic[ch])
	return dic

################### FUNC get_overlap_between_intervals  #####################################
#### Finds the overlap between two intervals end points inclusive
### a=[10,20]; b=[20,30] --> f(a,b)=1
### a=[10,20]; b=[15,30] --> f(a,b)=6
##############################################################################
def get_overlap_between_intervals(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0])+1)


