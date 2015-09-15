#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import gzip
import itertools
from itertools import islice
import numpy
import resource

REfilename=sys.argv[1] # /srv/gsfs0/projects/kundaje/users/mtaranov/projects/dynamic3D/FIT-HI-C/hic_flexible/data/reference_genomes/hg19/Digest_ENCODEHg19_HindIII_None_19-38-29_25-11-2014.txt
Resolution=int(sys.argv[2]) #10

outfilename=sys.argv[3] #test
outfilename2=sys.argv[4] #test

RE_dict={}
REFrag_dict={}
REboundary={}

# creates dictionary of individual RE sites: RE_dict={chr1: [[REstart1, REend1], [REstart2, REend2],...], chr2: [[]]}
for line in open(REfilename,'r'):
    words=line.rstrip().split()
    if words[0][:3]!="chr": # avoid header line
        continue
    if words[0][3:]=='M': #remove chrM
        continue
    if words[0] not in RE_dict:
        RE_dict[words[0]]=[[int(words[1]), int(words[2])]]
    else:
        RE_dict[words[0]].append([int(words[1]), int(words[2])])

# creates dictionary of RE boundaries: REboundary={chr1: [REend1, REend1, REend3,...], chr2:[REend1    , REend1, REend3,...], ...}
#    if words[0] not in REboundary:
#        REboundary[words[0]]=[int(words[2])]
#    else:
#        REboundary[words[0]].append(int(words[2]))


# creates dictionary fragments @RESOLUTION. For instance, at RE=10 {chr1: [[REstart1, REend10], [REstart11, REend20],...], chr2: [[]]} 
for key in RE_dict:
    NumberOfFrag=(len(RE_dict[key])/Resolution)*Resolution
    # if number of fragments in RE_dict can be devided by number of fragments @RESOLUTION without a reminder
    if  NumberOfFrag==len(RE_dict[key]):  
        for i in xrange(0, NumberOfFrag, Resolution):
            REFragCoor=[RE_dict[key][i][0], RE_dict[key][i+Resolution-1][1]]
            if key not in REFrag_dict:
                REFrag_dict[key]=[REFragCoor]
            else:
                REFrag_dict[key].append(REFragCoor)
    # if number of fragments in RE_dict can NOT be devided by number of fragments @RESOLUTION without a reminder, an extra RE fragment is created to sum up the remaining fragments
    else: 
        for i in xrange(0, NumberOfFrag+1, Resolution):
            if i==NumberOfFrag:
                REFragCoor=[RE_dict[key][i][0], RE_dict[key][len(RE_dict[key])-1][1]]
            else:
                REFragCoor=[RE_dict[key][i][0], RE_dict[key][i+Resolution-1][1]]
            if key not in REFrag_dict:
                REFrag_dict[key]=[REFragCoor]
            else:
               REFrag_dict[key].append(REFragCoor)
#    print key, len(REFrag_dict[key])

## writes out: chr# REstart  REend  RESiteMid  RESiteLength
outfile = open(outfilename,'w')
outfile2 = open(outfilename2,'w')
for key in REFrag_dict:
    for i in REFrag_dict[key]:
        outfile.write(key+'\t'+str(i[0])+'\t'+str(i[1])+'\t'+str(i[0]+(i[1]-i[0])/2)+'\t'+str(i[1]-i[0]+1)+'\n')
        outfile2.write(key+'\t'+str(i[0]+(i[1]-i[0])/2)+'\n')

