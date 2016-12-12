# written by Jin Lee, 2016

import os
import json
import subprocess
import collections

# find all qc_summary.json recursively
# json_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(os.getcwd()) \
#     for f in filenames if os.path.splitext(f)[1] == 'qc_summary.json']

json_files = subprocess.check_output("find . -name 'ENCODE_summary.json'", \
                                    shell=True ).strip().split('\n')
# read json
jsons = []
for json_file in json_files:
    # if not ('ENCSR765MXG' in json_file or 'ENCSR280ZDP' in json_file) : continue
    f = open(json_file,'r')
    jsons.append( json.load(f) )
    f.close()

with open('ENCODE_accession.csv','w') as f:
    
    # look at headers first
    raw_headers = list()

    for json in jsons:
        if not 'data_files' in json:
            continue
        data_files = json['data_files']
        for data_file in data_files:
            for key in data_file:
                if not key in raw_headers:
                    raw_headers.append( key )
    # sort header
    order_by_header = collections.defaultdict(int, \
        {
            'file_format':20,
            'file_format_type':19,
            'output_type':18,
            'dataset':17,
            'assembly':16,
            'aliases:array':15,
            'derived_from:array':14,
            'md5sum':13,
            'award':12,
            'lab':11,
            'submitted_file_name':10,
        })

    headers = sorted(raw_headers, key=lambda x: order_by_header[x], reverse=True)

    # write header
    f.write( ','.join( headers ) +'\n')

    lines = list()

    def find_submitted_file_name( submitted_file_name ):
        # recursively find file under a working directory and return path relative to working dir.
        files = subprocess.check_output("find . -type f -name '%s'" % (submitted_file_name), \
                    shell=True ).strip().split('\n')
        return files[0]

    # for each replicate, write contents
    for json in jsons:
        if not 'data_files' in json:
            continue
        data_files = json['data_files']
        for data_file in data_files:
            line = collections.OrderedDict()
            for key in headers:
                if key in data_file:
                    if key == 'submitted_file_name':
                        line[key] = find_submitted_file_name( data_file[key] )
                    else:
                        line[key] = data_file[key]
                else:
                    line[key] = ""
            lines.append(line)

    order_by_file_format = collections.defaultdict(int, \
        {
            'bam':20,
            'tagAlign':19,
            'bigWig':18,
            'bed':17,
            'bigBed':16,
        })
    order_by_output_type = collections.defaultdict(int, \
        {
            'alignments':20,
            'unfiltered alignments':19,
            'signal p-value':18,
            'fold change over control':17,
            'filtered peaks':16,
            'replicated peaks':15,
            'idr thresholded peaks':14,
            'optimal idr thresholded peaks':13,
            'conservative idr thresholded peaks':12,
        })
    # sort lines
    sorted_lines = sorted(lines, key = lambda x: (\
        order_by_file_format[x['file_format']],\
        order_by_output_type[x['output_type']]), reverse=True)

    for line in sorted_lines:
        result = ''
        for key in headers:
            result += (line[key]+ ('' if key==headers[-1] else ','))
        f.write( result + '\n' )
