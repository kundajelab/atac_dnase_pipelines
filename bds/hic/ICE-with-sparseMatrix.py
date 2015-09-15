#!/usr/bin/env python
'''
Created on Feb 18, 2013

@author: ferhat
'''

import sys
import numpy as np
import os

USAGE = """USAGE: ICE-with-sparseMatrix.py <contactCountsFile> <midsFile> <norm> <outfile> <lowMappThres>
Implements the ICE correction by Imakaev et al in sparse mode.
Can be slower than full matrix implementation but at least is feasible for bigger datasets

<contactCountsFile>	: file that includes the counts of interaction between locus pairs
		format : 		<chr ID 1>	<mid 1> 	<chr ID 2>	<mid 2>	<count>
						chr8	142250000	chr8	142350000	31

<midsFile>				: file with bin midpoints 
		format:
						<chr ID>	<mid>	<anythingElse> <mappabilityValue> <anythingElse>+ 
						chr10	50000	NA	0.65	...

<norm>					: l1 or l2 norm

<outfile>				: file to write the output to. This will be a file with a matrix of size <numberOfBins> X <numberOfBins> 

<lowMappThres>			: filter out loci that are less mappable than the threshold. Reasonable choice 0.5
"""
eps=1e-3
maxIter=1000
def main(contactCountsFile,midsFile,norm,outfilename,lowMappThres):

	(allFragsDic,allFragsDicReverse,badFrags)=assignIndicesToFrags(midsFile,lowMappThres)
	noOfFrags=0
	for ch in allFragsDic:
		noOfFrags+=len(allFragsDic[ch])
	print noOfFrags
	(Xvals,Xindices)=readContactsInSparseForm(contactCountsFile,allFragsDic,noOfFrags)
	print len(Xvals)
	#for i in range(len(Xvals)):
	#	print str(i)+"\t"+str(len(Xvals[i]))

	(XnewVals,XnewIndices,totalbias)=ICE(Xvals,Xindices,noOfFrags,norm,badFrags,allFragsDicReverse,outfilename+".biases")
	print len(Xvals)
	
	writeContactsInSparseForm(outfilename,allFragsDic,allFragsDicReverse,XnewVals,XnewIndices,badFrags)
	#writeContactsInSparseForm-notNecessary(contactCountsFile,outfilename,allFragsDic,allFragsDicReverse,XnewVals,XnewIndices,totalbias)
	#writeContactsInSparseForm-BUG-found-03-27-2014(outfilename,allFragsDic,allFragsDicReverse,XnewVals,XnewIndices)

	return

def	writeContactsInSparseForm(outfilename,allFragsDic,allFragsDicReverse,XnewVals,Xindices,badFrags):

	outfile=open(outfilename,'w')
	for indx1 in range(len(Xindices)):
		ch1,mid1=allFragsDicReverse[indx1]
		c=0
		for indx2 in Xindices[indx1]:
			if indx1>indx2 or indx1 in badFrags or indx2 in badFrags:# don't print twice
				c+=1
				continue
			ch2,mid2=allFragsDicReverse[indx2]
			newContactCount=float(XnewVals[indx1][c])
			c+=1
			outfile.write(str(ch1)+"\t"+str(mid1)+"\t"+str(ch2)+"\t"+str(mid2)+"\t"+str(newContactCount)+"\n")
	#
	outfile.close()
	return

#def	writeContactsInSparseForm-notNecessary(outfilename,allFragsDic,allFragsDicReverse,XnewVals,Xindices,totalbias):
#	infile=open(contactCountsFile,'r')
#	for line in infile:
#		ch1,mid1,ch2,mid2,contactCount=line.split()
#		contactCount=float(contactCount)
#		indx1=allFragsDic[ch1][mid1]
#		indx2=allFragsDic[ch2][mid2]
#		b1=totalbias[indx1]
#		b2=totalbias[indx2]
#		newContactCount=(1.0*contactCount)/(b1*b2)
#		outfile.write(str(ch1)+"\t"+str(mid1)+"\t"+str(ch2)+"\t"+str(mid2)+"\t"+str(newContactCount)+"\n")
#	#
#	outfile.close()
#	return

#def	writeContactsInSparseForm-BUG-found-03-27-2014(outfilename,allFragsDic,allFragsDicReverse,XnewVals,Xindices):
#	outfile=open(outfilename,'w')
#	for indx1 in range(len(Xindices)):
#		ch1,mid1=allFragsDicReverse[indx1]
#		c=0
#		for indx2 in Xindices[indx1]:
#			if indx1>indx2:# don't print twice
#				c+=1
#				continue
#			ch2,mid2=allFragsDicReverse[indx2]
#			newContactCount=float(XnewVals[indx1][c])
#			c+=1
#			outfile.write(str(ch1)+"\t"+str(mid1)+"\t"+str(ch2)+"\t"+str(mid2)+"\t"+str(newContactCount)+"\n")
#	#
#	outfile.close()
#	return

def ICE(Xvals,Xindices,noOfFrags,norm,badFrags,allFragsDicReverse,outfilename):
	XnewVals=Xvals
	XnewIndices=Xindices
	oldBias=None
	outfile=open(outfilename,'w')

	# filtering stage
	for i in range(noOfFrags): # This part is the reason for afterICE contact counts with BUG found-03-27-2014
		if i in badFrags:
			XnewVals[i]=[]
			XnewIndices[i]=[]
	#

	totalbias=1.0*np.ones((noOfFrags,))
	#np.reshape(np.array([1.0 for i in Xvals]),(noOfFrags, 1))
	for it in range(maxIter):
		if norm == 'l1':
			sumds = np.array([1.0*sum(i) for i in Xvals]).reshape((noOfFrags,))
		elif norm == 'l2':
			sumds = np.array([np.sqrt(sum(np.power(i,2))) for i in Xvals]).reshape((noOfFrags,))
		print np.shape(sumds)
		dbias=sumds
		print np.shape(dbias)
		print np.shape(totalbias)
		dbias /= dbias[dbias != 0].mean()
		dbias[dbias == 0] = 1
		print dbias[5:20]
		totalbias*=dbias
		print totalbias[5:20]
		# implement this line in sparse form --> X /= dbias.T * dbias
		for i in range(len(Xindices)):
			if i%10000==0 and i>0:
				print i
			biasI=dbias[i]
			c=0
			for j in Xindices[i]: # be careful this is not range(len(Xindices[i]))
				biasJ=dbias[j]
				XnewVals[i][c]= (1.0*XnewVals[i][c])/((1.0*biasI)*(1.0*biasJ))
				c+=1
			#
		#

		if oldBias is not None and np.abs(oldBias - dbias).sum() < eps:
			print "break at iteration %d" % (it,)
			break
		if oldBias is not None:
			print ('ICE at iteration %d %s' % (it, np.abs(oldBias - dbias).sum()))
		oldBias = dbias.copy()
	#
	for i in range(noOfFrags):
		ch,mid=allFragsDicReverse[i]
		outfile.write(str(ch) +"\t"+ str(mid)+"\t" +str(totalbias[i])+"\n")
	outfile.close()
	return (XnewVals,XnewIndices,totalbias)


def	readContactsInSparseForm(contactCountsFile,allFragsDic,noOfFrags):
	Xvals=[]
	Xindices=[]
	for i in range(noOfFrags):
		Xvals.append([])
		Xindices.append([])
	#
	infile=open(contactCountsFile,'r')
	for line in infile:
		ch1,mid1,ch2,mid2,contactCount=line.split()
		contactCount=float(contactCount)
		indx1=allFragsDic[ch1][mid1]
		indx2=allFragsDic[ch2][mid2]
		#print str(indx1)+"\t"+str(indx2)
		Xvals[indx1].append(contactCount)
		Xindices[indx1].append(indx2)
		Xvals[indx2].append(contactCount)
		Xindices[indx2].append(indx1)
	#
	infile.close()
	return (Xvals,Xindices)


def assignIndicesToFrags(midsFile,lowMappThres):
	badFrags=[]
	allFragsDic={}
	allFragsDicReverse={}
	infile=open(midsFile,'r')
	indx=0
	for line in infile:
		words=line.split()
		if words[0]=="chr": # avoid header line
			continue
		currChr=words[0]; currMid=words[1]; mapp=float(words[3]);
		if currChr not in allFragsDic:
			allFragsDic[currChr]={}
		allFragsDic[currChr][currMid]=indx
		allFragsDicReverse[indx]=[currChr,currMid]
		if mapp<=lowMappThres:
			badFrags.append(indx)
		indx+=1
	#END
	infile.close()
	return (allFragsDic,allFragsDicReverse,badFrags)

if __name__ == "__main__":
	if (len(sys.argv) != 6):
		sys.stderr.write(USAGE)
		sys.exit(1)
	main(str(sys.argv[1]),str(sys.argv[2]),str(sys.argv[3]),str(sys.argv[4]),float(sys.argv[5]))
