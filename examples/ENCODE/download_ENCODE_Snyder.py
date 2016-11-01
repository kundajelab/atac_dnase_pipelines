#!/usr/bin/env python2

import sys, json, urllib2

exps = ['ENCSR955JSO',\
'ENCSR355SGJ',\
'ENCSR846VLJ',\
'ENCSR386HAZ',\
'ENCSR404LLJ',\
'ENCSR761TKU',\
'ENCSR086OGH',\
'ENCSR668VCT',\
'ENCSR260ZIV',\
'ENCSR788TRR',\
'ENCSR670REK',\
'ENCSR337UIU',\
'ENCSR540BML',\
'ENCSR630REB',\
'ENCSR846ZBX',\
'ENCSR654UYP',\
'ENCSR078EBD',\
'ENCSR851SBY',\
'ENCSR548QCP']

fastq = dict()

for exp in exps:
    json_data = urllib2.urlopen('https://www.encodeproject.org/experiments/'+exp+'/?format=json').read()
    json_obj = json.loads(json_data)

    print '#============ exp: %s =========' % (exp,)
    for i,f in enumerate(json_obj['files']):
        out_dir = '$TMP/data/ENCODE_mouse/%s/pair%d' % (exp, i+1)
        print 'mkdir -p %s' % (out_dir,)
        print 'wget -P %s -bqc https://www.encodeproject.org%s' % ( out_dir, f['href'] )
