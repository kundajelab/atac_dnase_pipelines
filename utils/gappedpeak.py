#!/usr/bin/env python2

import sys,os

if len(sys.argv)!=3:
	print '<gappedpeak file> <track name>'
	sys.exit()

infile,outfile=sys.argv[1:]

id=1
fout=open(outfile,'w')
with open(infile) as fin:
	for line in fin:
		lst=line.rstrip().split('\t')
		fout.write('{0[0]}\t{0[1]}\t{0[2]}\tscorelst:[{0[12]},{0[13]},{0[14]}],id:{1},struct:{{thin:[[{0[1]},{0[2]}]],thick:['.format(lst,id))
		id+=1
		a=int(lst[1])
		sizes=lst[10].split(',')
		starts=lst[11].split(',')
		for i in range(len(sizes)):
			fout.write('[{0},{1}],'.format(a+int(starts[i]),a+int(starts[i])+int(sizes[i])))
		fout.write(']},')

		if len(lst[3])>1:
			fout.write('name:"'+lst[3]+'",')
		if lst[5]!='.':
			fout.write('strand:"'+lst[5]+'",')
		fout.write('\n')

fout.close()

os.system('sort -k1,1 -k2,2n '+outfile+' > '+outfile+'.srt')
os.system('mv '+outfile+'.srt'+' '+outfile)
os.system('bgzip -f '+outfile)
os.system('tabix -f -p bed '+outfile+'.gz')
