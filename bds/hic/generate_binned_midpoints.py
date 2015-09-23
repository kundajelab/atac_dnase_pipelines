#!/usr/bin/env python
'''
Created on Nov 30, 2011

@author: ferhat
'''

import sys
import math

USAGE = """USAGE: generate_binned_midpoints.py  <fragsize> <chrLengthsfilename> <outfilename>


"""


def generateBinnedMids(fragsize, chrLengthsfilename, outfilename):
    
    outFile=open(outfilename, 'w')

    chrLengthsFile=open(chrLengthsfilename, 'r')
    chrList=[]
    lengthsList=[]
    for line in chrLengthsFile:
        words=line.rstrip().split()
        #tempchr=(words[0])[3:]
        tempchr=(words[0])
        if tempchr!='M':
            #if tempchr=='X':
            #    tempchr='23'
            #elif tempchr=='Y':
            #    tempchr='24'
            chrList.append(tempchr)
            lengthsList.append(int(words[1]))
    chrLengthsFile.close()   

    for x in range(len(chrList)):
        chrLength = lengthsList[x]
        noOfFrags= math.ceil((chrLength-0.00001)/fragsize)
#        print noOfFrags
        for i in range(int(noOfFrags)):
            outFile.write(str(chrList[x])+'\t'+repr(fragsize/2+(i*fragsize)) +'\n')

    outFile.close() 
    
if __name__ == "__main__":
    if (len(sys.argv) != 4):
        sys.stderr.write(USAGE)
        sys.exit(1)
    generateBinnedMids(int(sys.argv[1]),str(sys.argv[2]),str(sys.argv[3]))
        
