#!/usr/bin/env python2.7

# Author: Jason Buenrostro, Stanford University
# The following program will compress daisy chain seq data into singe molecules

##### IMPORT MODULES #####
# import necessary for python
import os
import re
import sys
import gzip
import string
import Levenshtein
from optparse import OptionParser

##### DEFINE FUNCTIONS #####
# Reverse complement
complement = string.maketrans('ATCGN', 'TAGCN')
def reverse_complement(sequence):
    return sequence.upper().translate(complement)[::-1]

# Align with mismatch, find first and move on, assumes only one
def fuzz_align(s_seq,l_seq,mismatch):
    for i, base in enumerate(l_seq):  # loop through equal size windows
        l_subset = l_seq[i:i+len(s_seq)]
        dist = Levenshtein.distance(l_subset, s_seq)
        if dist <= mismatch:  # find first then break
            return i, dist
            break

# added by Jin Lee for hot fix (output name bug)
def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)

#### OPTIONS ####
# define options
opts = OptionParser()
usage = "usage: %prog [options] [inputs] This will trim adapters"
opts = OptionParser(usage=usage)
opts.add_option("-a", help="<Read1> Accepts fastq or fastq.gz")
opts.add_option("-b", help="<Read2> Accepts fastq or fastq.gz")
options, arguments = opts.parse_args()

# return usage information if no argvs given AND they're not available in the environment
# command line arguments always override environment variables
if len(sys.argv)==1:
	p1_in = os.environ.get('P1_IN')
	p2_in = os.environ.get('P2_IN')
	if  (p1_in is None) or (p2_in is None):
		os.system(sys.argv[0]+" --help")
		sys.exit()
else:
	##### INPUTS AND OUTPUTS #####
	# name input and outputs
	p1_in = options.a
	p2_in = options.b

# name outputs and print to working dir
p1_file = p1_in.split('/')[-1]
p2_file = p2_in.split('/')[-1]
p1_out = re.sub(".fastq", ".trim.fastq", p1_file)
p2_out = re.sub(".fastq", ".trim.fastq", p2_file)

#check for file type and open input file
append = p1_in.split('.')[-1]
if append == "fastq":
    p1_rds = open(p1_in,'r')
    p2_rds = open(p2_in,'r')
    p1_out = re.sub(".fastq", ".trim.fastq", p1_file)
    p2_out = re.sub(".fastq", ".trim.fastq", p2_file)
elif append == "fq":
    p1_rds = open(p1_in,'r')
    p2_rds = open(p2_in,'r')
    p1_out = re.sub(".fq", ".trim.fastq", p1_file)
    p2_out = re.sub(".fq", ".trim.fastq", p2_file)
elif append == "gz":
    p1_rds = gzip.open(p1_in,'r')
    p2_rds = gzip.open(p2_in,'r')
    p1_out = re.sub(".fastq.gz", ".trim.fastq", p1_file)
    p2_out = re.sub(".fastq.gz", ".trim.fastq", p2_file)
    p1_out = re.sub(".fq.gz", ".trim.fastq", p1_out)
    p2_out = re.sub(".fq.gz", ".trim.fastq", p2_out)
else:
    sys.exit("ERROR! The input file2 must be a .fastq or .fastq.gz")

##### SCRIPT #####
# initialize variables
i=0;j=0;k=0;tot_b=0;count=1
n=20  # match seq
mismatch=1  # only allow 0-1 mismatches for now, if allow two then gets mis indexed, to fix this need to change fuzz_align to save L as a vector and reiterate to find 2nd

# initilize write files
r1_write = open(p1_out, 'w')
r2_write = open(p2_out, 'w')

while 1:
    # read lines
    p1_line = p1_rds.readline()
    p2_line = p2_rds.readline()

    # break if at end of file
    if not p1_line:
        break

    # load fastq into memory
    if count ==1:
        seqhead1 = p1_line
        seqhead2 = p2_line
    elif count ==2:
        seq1 = p1_line.rstrip()
        seq2 = p2_line.rstrip()
    elif count ==3:
        qualhead1 = p1_line
        qualhead2 = p2_line
    elif count ==4:
        qual1 = p1_line.rstrip()
        qual2 = p2_line.rstrip()

        # align reads to themselves
        i = i+1  # total reads
        rc_seq2 = reverse_complement(seq2[0:n])
        idx = seq1.rfind(rc_seq2) # look for perfect match
        if idx > 0:
            j = j+1  # 0 mismatchs
        elif mismatch>0:
            hold = fuzz_align(rc_seq2,seq1,mismatch)  # else allow for mismatch
            if hold:
                idx,mis=hold
                if mis == 1:
                    k=k+1  # 1 mismatch

        # trim reads if idx exist
        if idx > 0:
            # keep track on how much trimming
            tot_b = tot_b+len(seq2[idx+n:-1]) #track total bases trimmed 
            
            # trim data
            seq1 = seq1[0:idx+n-1]  # modified to sub1 because some aligners (bowtie) dont like perfectly overlapping reads
            seq2 = seq2[0:idx+n-1]
            qual1 = qual1[0:idx+n-1]
            qual2 = qual2[0:idx+n-1]
        
        # print read1
        r1_write.write(seqhead1)
        r1_write.write(seq1+"\n")
        r1_write.write(qualhead1)
        r1_write.write(qual1+"\n")

        # print read2
        r2_write.write(seqhead2)
        r2_write.write(seq2+"\n")
        r2_write.write(qualhead2)
        r2_write.write(qual2+"\n")

    # increment count
    count = count + 1
    if count == 5:
        count = 1
    else:
        count = count

# close files to write the file
r1_write.close()
r2_write.close()
p1_rds.close()
p2_rds.close()

# write file output names for passing into next step of pipeline
# !!! DO NOT WRITE ANYTHING ELSE TO STDOUT AFTER THIS !!!
sys.stdout.write(p1_out + '\n')
sys.stdout.write(p2_out + '\n')

# give summary
sys.stderr.write(str(i)+" sequences total\n")
sys.stderr.write(str(j)+" sequences trimmed with 0 mismatches\n")
sys.stderr.write(str(k)+" sequences trimmed with 1 mismatch\n")
sys.stderr.write(str(tot_b/(j+k))+" mean number of bases trimmed for reads requiring trimming\n")
