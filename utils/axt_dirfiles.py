#!/usr/bin/env python2

import sys,glob,gzip,os

# axt format: http://genome.ucsc.edu/goldenPath/help/axt.html

if len(sys.argv)!=3:
	print '<chr size file> <output file> Run under the dir of gzipped Axt files, presumably one for each target chr but that doesn\'t matter'
	sys.exit()

chrsize={}
with open(sys.argv[1]) as fin:
	for line in fin:
		lst=line.rstrip().split('\t')
		chrsize[lst[0]]=int(lst[1])


OF=sys.argv[2]

fout=open(OF,'w')

id=1

for f in glob.glob('*'):
	fin=gzip.GzipFile(f,'r')
	line=fin.readline()
	while line:
		if line[0]!='#':
			lst=line.rstrip().split()
			# query start/stop
			a=0
			b=0
			if lst[7]=='+':
				a=int(lst[5])-1
				b=lst[6]
			else:
				c=chrsize[lst[4]]
				a=c-int(lst[6])
				b=c-int(lst[5])+1

			fout.write('{0[1]}\t{2}\t{0[3]}\tid:{1},genomealign:{{chr:"{0[4]}",start:{3},stop:{4},strand:"{0[7]}",targetseq:'.format(
				lst,
				id,
				int(lst[2])-1,
				a,
				b
				))
			id+=1
			line=fin.readline().rstrip()
			fout.write('"'+line+'",queryseq:')
			line=fin.readline().rstrip()
			fout.write('"'+line+'"}\n')
			fin.readline()
		line=fin.readline()


fout.close()


os.system('sort -k1,1 -k2,2n '+OF+' > xx')
os.system('mv xx '+OF)
os.system('bgzip -f '+OF)
os.system('tabix -f -p bed '+OF+'.gz')
