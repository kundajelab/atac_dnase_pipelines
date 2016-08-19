#!/usr/bin/env python2

# show -log10(GLOBAL IDR SCORE) instead of narrowpeak pval

import sys,os

if len(sys.argv)!=3:
	print '<narrowpeak_idr file> <track name>'
	sys.exit()

infile,outfile=sys.argv[1:]

id=1
fout=open(outfile,'w')
with open(infile) as fin:
	for line in fin:
		lst=line.rstrip().split('\t')
		fout.write('{0[0]}\t{0[1]}\t{0[2]}\tscorelst:[{0[6]},{0[7]},{0[8]},{0[10]},{0[11]}],id:{1},'.format(lst,id))
		id+=1
		if len(lst[3])>1:
			fout.write('name:"'+lst[3]+'",')
		else:
			fout.write('name:"'+str(id)+'",')
		if lst[5]!='.':
			fout.write('strand:"'+lst[5]+'",')
		if lst[9]!='-1':
			fout.write('sbstroke:['+lst[9]+']')
		fout.write('\n')

fout.close()

os.system('sort -k1,1 -k2,2n '+outfile+' > '+outfile+'.srt')
os.system('mv '+outfile+'.srt'+' '+outfile)
os.system('bgzip -f '+outfile)
os.system('tabix -f -p bed '+outfile+'.gz')
