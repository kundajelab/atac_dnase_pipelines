#!/usr/bin/env python

import sys
import gzip
#import pyparsing
#import matplotlib.pyplot as plt

window=sys.argv[1]
ContactPairFile=sys.argv[2]
BlacklistFile=sys.argv[3]
outfilename=sys.argv[4]
lib=sys.argv[5]

w=int(window)

BlackList=[]
# construct blacklist of windows 
for line in open(BlacklistFile, 'r'):
    chr=line.rstrip().split()[0]
    midpoint=line.rstrip().split()[1]
    BlackList.append([chr, midpoint])


Distance_black=[]
ListOfContacts_black=[]

outfile=open(outfilename,'w')
# pull only pairs which do not overlap Blacklist 
for line in open(ContactPairFile,'r'):
    words=line.rstrip().split()
    if ([words[0], words[1]] not in BlackList) and ([words[2], words[3]] not in BlackList):
        outfile.write(line)
#        # plot intra-chr interactions from the final list of contacts 
#        if words[0]==words[2]:
#            Distance_item_black=(int(words[3])-int(words[1]))/1000
#            Distance_black.append(Distance_item_black)
#            ListOfContacts_item_black=int(words[4])
#            ListOfContacts_black.append(ListOfContacts_item_black)   
#
#fig = plt.figure()
#plt.scatter(Distance_black, ListOfContacts_black, c='y')
#plt.xlim(40, 8000)
#plt.ylim(0, 2000)
#plt.xlabel('Genomic Distance (kbp)')
#plt.ylabel('Contact Counts')
#fig.suptitle('Intra-chr Contact Counts outside the Blacklist '+str(lib)+' w='+str(w))
#fig.savefig(outfilename)     
