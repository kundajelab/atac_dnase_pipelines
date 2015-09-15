from optparse import OptionParser
import os
import math
from time import gmtime, strftime
import re
import gzip
import random
import numpy 
'''
Author:Oana Ursu
'''

def main():
    parser=OptionParser()
    
    parser.add_option('--read_length',dest='readLen',default='50',help='Length of the reads you are using.')
    parser.add_option('--map_sufix',dest='map_suf',default='.uint8.unique',help='Default is set for Kundajelab to ".uint8.unique"')
    parser.add_option('--map_prefix',dest='map_prefix',default='/srv/gsfs0/projects/kundaje/commonRepository/mappability/encodeHg19Male/globalmap_k20tok54/chr',help='Default is already set for KUndajelab to be "/srv/gsfs0/projects/kundaje/commonRepository/mappability/encodeHg19Male/globalmap_k20tok54/chr"')
    parser.add_option('--regions',dest='bed',help='Bed file with regions for which we want mappability, independent of strand. Only the first 3 entries matter. The rest of the entries will be copied over from the "bed file"',default='/home/oursu/testbed')
    parser.add_option('--out',dest='out',help='Outfile',default='/home/oursu/maptest')
    opts,args=parser.parse_args()

    #output file
    out=open(opts.out,'w')
    readLength=int(opts.readLen)

    cur_chromo=''
    for line in open(opts.bed,'r').readlines():
        items=line.strip().split('\t')
        chromo=re.sub('chr','',items[0])
        start=int(items[1])
        end=int(items[2])
        other=items[3:]
        mid=int(float((start+end)/2))
        if chromo!=cur_chromo:
            fname=opts.map_prefix+chromo+opts.map_suf
            print fname
            cur_chromo=chromo
            #mappability on plus strand
            mappability_vector=numpy.fromfile(open(fname,'rb'), dtype=numpy.uint8)
            notZero=(mappability_vector!=0)*1
            uniquelyMappableBelowThreshold=(mappability_vector<=readLength)*1
            #mappability on minus strand
            negStrandMappability=0*mappability_vector
            negStrandMappability[readLength:len(negStrandMappability)]=mappability_vector[0:(len(negStrandMappability)-readLength)]
            notZeroNeg=(negStrandMappability!=0)*1
            uniquelyMappableBelowThresholdNeg=(negStrandMappability<=readLength)*1
            #put everything together in mappability scores - 1 if mappable for plus and minus strand
            mappability_scores=numpy.multiply(numpy.multiply(notZero,uniquelyMappableBelowThreshold),numpy.multiply(notZeroNeg,uniquelyMappableBelowThresholdNeg))
        practical_end=int(min(end,len(mappability_scores)))
        cur_mappability=float(sum(mappability_scores[start:practical_end]))/float(abs((start-practical_end)))
        out.write('chr'+chromo+'\t'+str(mid)+'\t'+'1'+'\t'+str(cur_mappability)+'\n')
    out.close()

main()
