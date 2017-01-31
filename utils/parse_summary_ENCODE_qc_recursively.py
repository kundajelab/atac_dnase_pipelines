#!/usr/bin/env python2

# written by Jin Lee, 2016

import os, sys
import json
import subprocess
import collections
import argparse

parser = argparse.ArgumentParser(prog='ENCODE_summary.json parser for ENCODE QC import', \
                                    description='Recursively find ENCODE_summary.json, parse it and make an excel file for importing quality metrics to the ENCODE portal. Use https://github.com/ENCODE-DCC/pyencoded-tools/blob/master/ENCODE_import_data.py for uploading.')
parser.add_argument('out-file-prefix', type=str, \
                        help='Output TSV filename prefix, output files will be [PREFIX].[QC_TYPE].tsv')
parser.add_argument('--search-dir', type=str, default='.', \
                        help='Root directory to search for ENCODE_summary.json')
parser.add_argument('--json-file', type=str, default='ENCODE_summary.json', \
                        help='Specify json file name to be parsed')
parser.add_argument('--sort-by-genome-and-exp', dest='sort_by_genome_and_exp', action='store_true', \
                        help='Sort rows by genomes and ENCODE experiment accession ID')
parser.add_argument('--ignored-accession-ids-file', type=str, \
                        help='Accession IDs in this text file will be ignored. (1 acc. ID per line)')
parser.set_defaults(sort_by_genome_and_exp=False)

args = parser.parse_args()

# loaded ignored accession list
ignored_accession_ids = []
if args.ignored_accession_ids_file and os.path.isfile(args.ignored_accession_ids_file):
    with open(args.ignored_accession_ids_file,'r') as f:
        ignored_accession_ids = f.read().splitlines()
    ignored_accession_ids = \
        [accession_id for accession_id in ignored_accession_ids if accession_id and not accession_id.startswith("#") ]

# find all ENCODE_summary.json recursively
json_files = subprocess.check_output("find %s -name %s" % (args.search_dir,args.json_file), \
                                    shell=True ).strip().split('\n')
# read json
jsons = []
for json_file in json_files:
    with open(json_file,'r') as f:
        jsons.append( json.load(f) )
 
# look at headers first
raw_headers = list()

for json in jsons:
    if not 'ENCODE_quality_metrics' in json:
        continue
    data_files = json['ENCODE_quality_metrics']
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
args.out_file_prefix.write( ','.join( headers ) +'\n')

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
    if json['ENCODE_accession'] in ignored_accession_ids: continue
    data_files = json['ENCODE_quality_metrics']
    for data_file in data_files:
        line = collections.OrderedDict()
        for key in headers:
            if key in data_file:
                # if key == 'submitted_file_name':
                #     line[key] = find_submitted_file_name( data_file[key] )
                # else:
                line[key] = data_file[key]
            else:
                line[key] = ""
        lines.append(line)

# sort lines
if args.sort_by_genome_and_exp:
    sorted_lines = sorted(lines, key = lambda x: (\
        x['assembly'],\
        x['dataset']) )
else:
    sorted_lines = lines

for line in sorted_lines:
    result = ''
    for key in headers:
        result += (line[key]+ ('' if key==headers[-1] else ','))
    args.out_file_prefix .write( result + '\n' )
