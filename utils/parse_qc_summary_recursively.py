# written by Jin Lee, 2016

import os
import json
import subprocess

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
sorted_jsons = sorted(jsons, key = lambda x: (\
    x['ENCODE_award_rfa'], \
    x['ENCODE_assay_category'], \
    x['ENCODE_assay_title'], \
    x['species'], \
    x['title']))

f = open('qc_summary.tsv','w')
    
# look at headers first
headers = [( 'common',\
        'ENCODE award rfa\t'+\
        'ENCODE assay category\t'+\
        'ENCODE assay title\t'+\
        'species\t'+\
        'title\t'+\
        'replicate')]

for json in jsons:
    for qc_file in json['qc_files']:
        pair = ( qc_file['qc_type'], qc_file['header'] )
        if not pair in headers:
            headers.append( pair )

# write header1
f.write( '\t'.join( [ qc_type+'\t'*(len(header_line.split('\t'))-1) \
                        for qc_type, header_line in headers ] ) +'\n')

# write header2
f.write( '\t'.join( [ header_line \
                        for qc_type, header_line in headers ] ) +'\n')

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
        for qc_type, header_line in headers:
            if qc_type=='common':
                continue
            found = False
            for qc_file in json['qc_files']:
                if rep != qc_file['info']:
                    continue
                if qc_type == qc_file['qc_type']:
                    result += ('\t'+qc_file['contents'] )
                    found = True
                    break
            if not found:
                result+= ('\t'+'\t'*len(header_line.split('\t')) )

        f.write( result + '\n' )

f.close()