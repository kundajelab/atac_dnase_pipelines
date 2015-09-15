#!/usr/bin/env python
import sys
import gzip
import matplotlib as mpl
import matplotlib.pyplot as plt
import math

BiaseFile=sys.argv[1]
outfilename=sys.argv[2]
lib=sys.argv[3]
w=sys.argv[4]

#data = [math.log1p(float(line.rstrip().split()[2])) for line in gzip.open(BiaseFile)]
data = [float(line.rstrip().split()[2]) for line in gzip.open(BiaseFile)]

fig1 = plt.figure(1)
plt.hist(data)
plt.xlabel('Bias Value')
fig1.suptitle(str(lib)+' w='+str(w))
plt.savefig(outfilename)
