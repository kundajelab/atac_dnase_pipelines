#!/usr/bin/env python2
# code extracted from Daniel Kim's ATAQC module (run_ataqc.py)

import os, sys, re, gzip

def getFileHandle(filename, mode="r"):
    if (re.search('.gz$',filename) or re.search('.gzip',filename)):
        if (mode=="r"):
            mode="rb";
        return gzip.open(filename,mode)
    else:
        return open(filename,mode)

def get_read_length(fastq_file):
    '''
    Get read length out of fastq file
    '''
    total_reads_to_consider = 1000000
    line_num = 0
    total_reads_considered = 0
    max_length = 0
    with getFileHandle(fastq_file, 'rb') as fp:
        for line in fp:
            if line_num % 4 == 1:
                if len(line.strip()) > max_length:
                    max_length = len(line.strip())
                total_reads_considered += 1
            if total_reads_considered >= total_reads_to_consider:
                break
            line_num += 1

    return int(max_length)

def main():
    print get_read_length(sys.argv[1])
    
if __name__ == "__main__":
    main()
