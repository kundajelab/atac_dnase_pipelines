#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import gzip

infilename=sys.argv[1]
lib=sys.argv[2]
w=sys.argv[3]
outfilename=sys.argv[4]

words=[line.rstrip().split() for line in gzip.open(infilename,'r')]

outfile=open(outfilename,'a')
Distance=[]
ListOfContacts=[]
Distance_vs_ListOfContacts=[]
for i in words:
    if i[0]==i[2]:
        Distance_item=(int(i[3])-int(i[1]))/1000
        Distance.append(Distance_item)
        ListOfContacts_item=int(i[4])
        ListOfContacts.append(ListOfContacts_item)
        outfile.write(str(Distance_item)+'\t'+str(ListOfContacts_item)+'\n')

fig = plt.figure()
plt.scatter(Distance, ListOfContacts)
plt.xlim(40, 8000)
plt.ylim(0, 2000)
fig.suptitle('Contact Counts '+str(lib)+' at RE'+str(w))
plt.xlabel('Genomic Distance (kbp)')
plt.ylabel('Contact Counts')
fig.savefig(outfilename)
