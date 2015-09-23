#!/usr/bin/env python

import sys
import gzip
import itertools

#mapThres=0.5
window=sys.argv[1]
BlacklistFile=sys.argv[2]
outfilename=sys.argv[3]

w=int(window)
BlackList=[] # blacklist of windows which overlap with blacklist motifs
# construct blacklist of windows 
for line in gzip.open(BlacklistFile, 'r'): # for each motif in blacklist
    chr=line.rstrip().split()[0]
    MotifStart=int(line.rstrip().split()[1])
    MotifEnd=int(line.rstrip().split()[2])
    MotifLength=(MotifEnd-MotifStart)/w
    # midpoint of the window that overlaps with the start of the blacklist motif
    StartMid=(MotifStart/w)*w+w/2
    # append window (window's midpoint) which overlaps with the start of the blacklist motif  
    BlackList.append([chr, str(StartMid)])
    # midpoint of the window that overlaps with the end of the blacklist motif
    EndMid=(MotifEnd/w)*w+w/2
    # number of windows that overlap with the current motif
    NumOfWinds_btw_StartMid_and_EndMid=(EndMid-StartMid)/w
    n=1
    # append windows (window's midpoints) which overlaps with the current blacklist motif
    while n<=NumOfWinds_btw_StartMid_and_EndMid:
        BlackList.append([chr, str(StartMid+w*n)])
        n = n + 1
# pulls unique windows 
BlackList.sort()
BlackList=list(BlackList for BlackList,_ in itertools.groupby(BlackList))

outfile=open(outfilename,'w')
for row in BlackList:
    outfile.write(" ".join(map(str,row))+"\n")
