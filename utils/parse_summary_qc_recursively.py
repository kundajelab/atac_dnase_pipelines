#!/usr/bin/env python2

# written by Jin Lee, 2016

import os
import sys
import re
import argparse
import json
import subprocess
from collections import OrderedDict

parser = argparse.ArgumentParser(prog='ENCODE_summary.json parser for QC', \
                                    description='Recursively find ENCODE_summary.json, parse it and make a TSV spreadsheet of QC metrics.')
parser.add_argument('--out-file', type=argparse.FileType('w'), default=sys.stdout, \
                        help='Output TSV filename)')
parser.add_argument('--search-dir', type=str, default='.', \
                        help='Root directory to search for ENCODE_summary.json')
parser.add_argument('--json-file', type=str, default='ENCODE_summary.json', \
                        help='Specify json file name to be parsed')

args = parser.parse_args()

# find all qc_summary.json recursively
# json_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(os.getcwd()) \
#     for f in filenames if os.path.splitext(f)[1] == 'qc_summary.json']

# find all ENCODE_summary.json recursively
json_files = subprocess.check_output("find -L %s -name %s" % (args.search_dir,args.json_file), \
                                    shell=True ).strip().split('\n')
# read json
jsons = []
for json_file in json_files:
    with open(json_file,'r') as f:
        jsons.append( json.load(f, object_pairs_hook=OrderedDict) )

# sort
# sorted_jsons = sorted(jsons, key = lambda x: (\
#     x['ENCODE_award_rfa'], \
#     x['ENCODE_assay_category'], \
#     x['ENCODE_assay_title'], \
#     x['species'], \
#     x['title']))

# look at headers first
headers = OrderedDict()
headers['common'] = [\
        'ENCODE award rfa',\
        'ENCODE assay category',\
        'ENCODE assay title',\
        'species',\
        'title',\
        'replicate']

# first take longest header for each qc_type
for json in jsons:
    for qc_file in json['qc_files']:
        qc_type = qc_file['qc_type']
        if qc_type == 'pbc_PE':
            qc_type = 'pbc'
            qc_file['qc_type'] = qc_type
        header_list = qc_file['header'].split('\t')        
        if not qc_type in headers or len(headers[qc_type])<len(header_list):
            headers[qc_type] = header_list

qc_type = 'files_to_be_submitted'
headers[qc_type] = []

# second add missing items for each qc_type
for json in jsons:
    for qc_file in json['qc_files']:
        qc_type = qc_file['qc_type']
        header_list = qc_file['header'].split('\t')
        for header_item in header_list:
            if not header_item in headers[qc_type]:
                headers[qc_type].append(header_item)
    # files to be submitted to ENCODE portal
    qc_type = 'files_to_be_submitted'
    for data_file in json['data_files']:
        header_item = ":".join([data_file['output_type'],data_file['file_format']])
	if not header_item in headers[qc_type]:
	    headers[qc_type].append(header_item)

# write header1
args.out_file.write( '\t'.join( [ qc_type+'\t'*(len(headers[qc_type])-1) \
                        for qc_type in headers ] ) +'\n')

# write header2
headers_wo_numbering = OrderedDict()
for qc_type in headers:
    headers_wo_numbering[qc_type] = [re.sub(r'^\d+_','',header) for header in headers[qc_type]]
args.out_file.write( '\t'.join( [ '\t'.join(headers_wo_numbering[qc_type]) \
                        for qc_type in headers_wo_numbering ] ) +'\n')

# for each replicate, write contents
for json in jsons:
    # count # of replicates per sample
    replicates = set()
    for qc_file in json['qc_files']:        
        info = qc_file['info'].replace('-pr','' )
        if not info or info == 'null': info = 'rep1'
        if not re.match(r'^rep\d+$', info): continue
        replicates.add( info )

    for rep in sorted(replicates):
        result = json['ENCODE_award_rfa']+'\t'+\
            json['ENCODE_assay_category']+'\t'+\
            json['ENCODE_assay_title']+'\t'+\
            json['species']+'\t'+\
            json['title']+'\t'+\
            rep
        for qc_type in headers:
    	    if rep == 'rep1' and qc_type == 'files_to_be_submitted':
		for header in headers[qc_type]:
                    header_found = False
		    tmp_result = []
    		    for data_file in json['data_files']:
                        if header == ':'.join([data_file['output_type'],data_file['file_format']]):
                            tmp_result.append(data_file['submitted_file_name'])
                            header_found = True
		    if header_found:
                        result += '\t'+ ','.join(tmp_result)		
                    else:
                        result += '\t'

            if qc_type=='common':
                continue
            registered_header_list = headers[qc_type]
            found = False
            for qc_file in json['qc_files']:
                info = qc_file['info'].replace('-pr','' )
                if not info or info == 'null': info = 'rep1'
                if not re.match(r'^rep\d+$', info): continue
                if rep != info:
                    continue
                if qc_type == qc_file['qc_type']:
                    header_list = qc_file['header'].split('\t')
                    contents_list = qc_file['contents'].split('\t')
                    h_to_c = dict(zip(header_list,contents_list))
                    for header_item in registered_header_list:
                        if header_item in header_list:
                            result += ('\t'+h_to_c[header_item] )
                        else:
                            result += ('\t')
                    found = True
                    break
		
            if not found:
                result += ('\t'*len(registered_header_list))

        args.out_file.write( result + '\n' )
