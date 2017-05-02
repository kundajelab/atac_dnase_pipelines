#!/usr/bin/env python2

# written by Jin Lee, 2016

import os, sys
import json
import subprocess
import collections
import argparse

parser = argparse.ArgumentParser(prog='ENCODE_summary.json parser for ENCODE accession', \
                                    description='Recursively find ENCODE_summary.json, parse it and make a CSV for uploading to the ENCODE portal. Use https://github.com/ENCODE-DCC/pyencoded-tools/blob/master/ENCODE_submit_files.py for uploading.')
parser.add_argument('--out-file', type=argparse.FileType('w'), default=sys.stdout, \
                        help='Output CSV filename)')
parser.add_argument('--search-dir', type=str, default='.', \
                        help='Root directory to search for ENCODE_summary.json')
parser.add_argument('--json-file', type=str, default='ENCODE_summary.json', \
                        help='Specify json file name to be parsed')
parser.add_argument('--sort-by-genome-and-exp', dest='sort_by_genome_and_exp', action='store_true', \
                        help='Sort rows by genomes and ENCODE experiment accession ID')
group_accession_ids = parser.add_mutually_exclusive_group()
group_accession_ids.add_argument('--ignored-accession-ids-file', type=str, \
                        help='Accession IDs in this text file will be ignored. (1 acc. ID per line)')
group_accession_ids.add_argument('--accession-ids-file', type=str, \
                        help='Only accession IDs in this text file will be downloaded. (1 acc. ID per line). Others will be ignored.')
parser.set_defaults(sort_by_genome_and_exp=False)

args = parser.parse_args()

# loaded ignored accession list
ignored_accession_ids = []
if args.ignored_accession_ids_file and os.path.isfile(args.ignored_accession_ids_file):
    with open(args.ignored_accession_ids_file,'r') as f:
        ignored_accession_ids = f.read().splitlines()
    ignored_accession_ids = \
        [accession_id for accession_id in ignored_accession_ids if accession_id and not accession_id.startswith("#") ]
accession_ids = []
if args.accession_ids_file and os.path.isfile(args.accession_ids_file):
    with open(args.accession_ids_file,'r') as f:
        accession_ids = f.read().splitlines()
    accession_ids = \
        [accession_id for accession_id in accession_ids if accession_id and not accession_id.startswith("#") ]

# find all ENCODE_summary.json recursively
json_files = subprocess.check_output("find -L %s -name %s" % (args.search_dir,args.json_file), \
                                    shell=True ).strip().split('\n')
# read json
jsons = []
for json_file in json_files:
    with open(json_file,'r') as f:
        jsons.append( json.load(f) )
 
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
args.out_file.write( ','.join( headers ) +'\n')

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
    if ignored_accession_ids and json['ENCODE_accession'] in ignored_accession_ids: continue
    if accession_ids and not json['ENCODE_accession'] in accession_ids: continue
    data_files = json['data_files']
    for data_file in data_files:
        line = collections.OrderedDict()
        for key in headers:
            if key in data_file:
                if key == 'submitted_file_name':
                    line[key] = find_submitted_file_name( data_file[key] )
                    metadata_file = line[key]+'.meta'
                    if os.path.exists(metadata_file):
                        with open(metadata_file,mode='r') as f:
                            for l in f:
                                if 'md5sum' in l:
                                    md5sum = l.split('=')[1].strip()
                                    if md5sum != data_file['md5sum']:
                                        print('Warning: In accession {}, md5sum of file {} does not match! (json:{}, actual:{})'.format(
                                            json['ENCODE_accession'], line[key], data_file['md5sum'], md5sum ))
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

if args.sort_by_genome_and_exp:
    sorted_lines = sorted(sorted_lines, key = lambda x: (\
        x['assembly'],\
        x['dataset']) )

for line in sorted_lines:
    result = ''
    for key in headers:
        result += (line[key]+ ('' if key==headers[-1] else ','))
    args.out_file.write( result + '\n' )
