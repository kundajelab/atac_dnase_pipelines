#!/usr/bin/env python2

# written by Jin Lee, 2016

import os, sys
import json
import subprocess
import collections
import argparse
import xlwt

parser = argparse.ArgumentParser(prog='ENCODE_summary.json parser for ENCODE QC import', \
                                    description='Recursively find ENCODE_summary.json, parse it and make an excel file for importing quality metrics to the ENCODE portal. Use https://github.com/ENCODE-DCC/pyencoded-tools/blob/master/ENCODE_import_data.py for uploading.')
parser.add_argument('out_file', metavar='out-file', type=str, \
                        help='Output Excel filename (extention should be .xls, not .xlsx)')
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
                        help='Only accession IDs in this text file will be parsed. (1 acc. ID per line). Others will be ignored.')
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
raw_headers = dict()

for json in jsons:
    if ignored_accession_ids and json['ENCODE_accession'] in ignored_accession_ids: continue
    if accession_ids and not json['ENCODE_accession'] in accession_ids: continue

    if not 'ENCODE_quality_metrics' in json: continue
    data_files = json['ENCODE_quality_metrics']
    for data_file in data_files:
        print data_file
        ENCODE_qc_type = data_file["ENCODE_qc_type"]
        if not raw_headers.has_key( "ENCODE_qc_type" ):
            raw_headers[ ENCODE_qc_type ] = list()
        for key in data_file:
            if key == "ENCODE_qc_type": continue
            if not key in raw_headers[ ENCODE_qc_type ]:
                raw_headers[ ENCODE_qc_type ].append( key )

# write header (fhs=file handles)
workbook = xlwt.Workbook()
sheets = {}

cnt=0
for ENCODE_qc_type in raw_headers:
    title = "".join([word.title().replace("Idr","IDR") for word in ENCODE_qc_type.split("_")])
    print "Creating a sheet with name: ", title
    # sheet = workbook.add_sheet(str(cnt))
    sheet = workbook.add_sheet(title)
    sheets[ENCODE_qc_type] = sheet
    for i, header in enumerate(raw_headers[ENCODE_qc_type]):
        sheet.write(0,i,header)
    cnt+=1
    # fh = open( "%s.%s.tsv" % (args.out_file_prefix,ENCODE_qc_type) ,'w')
    # fh.write(delimiter.join(raw_headers[ENCODE_qc_type]))
    # fh.write("\n")
    # fhs[ENCODE_qc_type] = fh

# for each replicate, write contents
lines = dict()
for json in jsons:
    if ignored_accession_ids and json['ENCODE_accession'] in ignored_accession_ids: continue
    if accession_ids and not json['ENCODE_accession'] in accession_ids: continue
    data_files = json['ENCODE_quality_metrics']
    for data_file in data_files:
        ENCODE_qc_type = data_file["ENCODE_qc_type"]
        if not lines.has_key(ENCODE_qc_type): 
            lines[ENCODE_qc_type] = list()
        line = collections.OrderedDict()
        for key in raw_headers[ENCODE_qc_type]:
            if key in data_file:
                line[key] = data_file[key]
            else:
                line[key] = ""
        lines[ENCODE_qc_type].append(line)

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def is_bool(s):
    if s.lower() in ['true','t','false','f']:
        return True
    else:
        return False

sorted_lines = lines
for ENCODE_qc_type in sorted_lines:
    data = sorted_lines[ENCODE_qc_type]
    sheet = sheets[ENCODE_qc_type]
    row = 1
    for line in data:
        for col, key in enumerate(line):
            val = line[key]
            if key.endswith('_pct'):
                val += "%"
            if val.startswith('null') or val.startswith('N/A:N/A'):
                val = "null"
            if is_int(val):
                val = int(val)
            elif is_float(val):
                val = float(val)
            # elif is_bool(val):
            # else:                
            #     style = xlwt.easyxf()
            sheet.write(row, col, label=val)
        row += 1

workbook.save(args.out_file)
