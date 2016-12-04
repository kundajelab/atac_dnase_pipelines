# written by Jin Lee, 2016

import os
import json
import subprocess
import collections

# find all qc_summary.json recursively
# json_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(os.getcwd()) \
#     for f in filenames if os.path.splitext(f)[1] == 'qc_summary.json']

json_files = subprocess.check_output("find . -name 'qc_summary.json'", \
                                    shell=True ).strip().split('\n')
# read json
jsons = []
for json_file in json_files:
    f = open(json_file,'r')
    jsons.append( json.load(f) )
    f.close()

# sort
# sorted_jsons = sorted(jsons, key = lambda x: (\
#     x['ENCODE_award_rfa'], \
#     x['ENCODE_assay_category'], \
#     x['ENCODE_assay_title'], \
#     x['species'], \
#     x['title']))

f = open('qc_summary.tsv','w')
    
# look at headers first
headers = collections.OrderedDict()
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
        header_list = qc_file['header'].split('\t')        
        if not qc_type in headers or len(headers[qc_type])<len(header_list):
            headers[qc_type] = header_list

# second add missing items for each qc_type
for json in jsons:
    for qc_file in json['qc_files']:
        qc_type = qc_file['qc_type']
        header_list = qc_file['header'].split('\t')
        for header_item in header_list:
            if not header_item in headers[qc_type]:
                headers[qc_type].append(header_item)

# write header1
f.write( '\t'.join( [ qc_type+'\t'*(len(headers[qc_type])-1) \
                        for qc_type in headers ] ) +'\n')

# write header2
f.write( '\t'.join( [ '\t'.join(headers[qc_type]) \
                        for qc_type in headers ] ) +'\n')

# for each replicate, write contents
for json in jsons:
    # count # of replicates per sample
    replicates = set()
    for qc_file in json['qc_files']:
        replicates.add( qc_file['info'] )

    for rep in sorted(replicates):
        result = json['ENCODE_award_rfa']+'\t'+\
            json['ENCODE_assay_category']+'\t'+\
            json['ENCODE_assay_title']+'\t'+\
            json['species']+'\t'+\
            json['title']+'\t'+\
            rep
        for qc_type in headers:
            if qc_type=='common':
                continue
            registered_header_list = headers[qc_type]
            found = False
            for qc_file in json['qc_files']:
                if rep != qc_file['info']:
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

        f.write( result + '\n' )

f.close()