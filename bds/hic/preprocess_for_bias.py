from optparse import OptionParser
import os
import math
from time import gmtime, strftime
import re
import gzip
import random
import numpy 
import sys
from Bio import SeqIO 
'''
Author:Oana Ursu
'''

def main():
    parser=OptionParser()
    parser.add_option('--chrSizes',dest='chrSizes',default='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/ref.fa.fai')
    parser.add_option('--genomedir',dest='gdir',default='/srv/gs1/projects/kundaje/oursu/Alignment/data/ENCODE_genomes/male/')
    parser.add_option('--REsite_file',dest='REfile',default='/srv/gsfs0/projects/kundaje/users/oursu/LongRangeInteractionPrediction/data/HICUP_genomes/Digest_ENCODEHg19_HindIII_None_19-38-29_25-11-2014.txt',help='A file with RE sites, as outputted from HICUP')
    parser.add_option('--read_length',dest='readLen',default='50',help='Length of the reads you are using.')
    parser.add_option('--map_sufix',dest='map_suf',default='.uint8.unique',help='Default is set for Kundajelab to ".uint8.unique"')
    parser.add_option('--map_prefix',dest='map_prefix',default='/srv/gsfs0/projects/kundaje/commonRepository/mappability/encodeHg19Male/globalmap_k20tok54/chr',help='Default is already set for KUndajelab to be "/srv/gsfs0/projects/kundaje/commonRepository/mappability/encodeHg19Male/globalmap_k20tok54/chr"')
    parser.add_option('--fragments',dest='frags',help='How to combine frags. Can specify a fixed size window in kb (e.g. 10kb) or a size in RE units (e.g. 10RE)')
    parser.add_option('--out',dest='out',help='Outfile',default='/home/oursu/maptest')
    opts,args=parser.parse_args()

    #output file
    out=open(opts.out,'w')

    #go through the RE file, and in each chromosome, combine fragments. remember them in a dictionary.
    regions={}
    #combine multiple fragments into one
    if 'RE' in opts.frags:
        numFrags=int(re.sub('RE','',opts.frags))
        print numFrags
        for line in open(opts.REfile,'r').readlines()[1:20]:
            print line
            items=line.strip().split('\t')
            chromo=re.sub('chr','',items[0])
            start=items[1]
            end=items[2]
            if chromo not in regions.keys():
                regions[chromo]={}
                fragCounter=0
                start_coord=start
            if fragCounter==0:
                start_coord=start
            fragCounter=fragCounter+1
            regions[chromo][start_coord]={}
            regions[chromo][start_coord]['end']=end
            regions[chromo][start_coord]['NumREsites']=fragCounter
            if fragCounter==numFrags:
                fragCounter=0
        print regions
    
    if 'kb' in opts.frags:
        wsize=int(1000*float(re.sub('kb','',opts.frags)))
        print wsize
        for line in open(opts.chrSizes,'r').readlines():
            print line
            items=line.strip().split('\t')
            chromo=items[0]
            chrsize=items[1]
            chr_pos=0
            regions[chromo]={}
            while chr_pos<chrsize:
                regions[chromo][str(chr_pos+1)]=chr_pos+wsize
                chr_pos=chr_pos+wsize
            print regions

    for chromo in regions.keys():
        try:
            #open chromosome mappability vector
            fname=opts.map_prefix+chromo+opts.map_suf
            mappability_vector=numpy.fromfile(open(fname,'rb'), dtype=numpy.uint8)
            notZero=(mappability_vector!=0)*1
            uniquelyMappableBelowThreshold=(mappability_vector<=opts.readLen)*1
            mappability_scores=numpy.multiply(notZero,uniquelyMappableBelowThreshold)
            #get all regions for this chromosome, and assign them mappabilities
            #for gc fafile=opts.gdir+'chr'+chromo+'.fa'
            #for gc faseq=SeqIO.parse(fafile, "fasta") 
            for region in regions[chromo].keys():
                start=int(region)
                end=int(regions[chromo][region]['end'])
                cur_mappability=float(sum(mappability_scores[start:end]))/float(abs((start-end)))
                regions[chromo][region]['map']=cur_mappability
                #for gc sequence=faseq[start:end]
                #for gc total_GC = sequence.count('G') + sequence.count('C')
                #for gc regions[chromo][region]['gc']=(float(total_GC) / float(abs(end-start)))
        except:
            print 'error on chromosome '+chromo
    print regions

    #write out including midpoint (all i still need is gc)
    for chromo in regions.keys():
        if chromo=='Chromosome':
            continue
        for region in regions[chromo].keys():
            midpoint=int(float(int(region)+int(regions[chromo][region]['end']))/2)
            out.write(chromo+'\t'+str(midpoint)+'\t'+str(regions[chromo][region]['NumREsites'])+'\t'+str(regions[chromo][region]['map'])+'\n') 

    out.close()

main()
