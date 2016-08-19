#!/usr/bin/env python2

import sys

if len(sys.argv)!=3:
	print '<input coordinate file> <basename of output file>'
	sys.exit()

infile,outn=sys.argv[1:]

aliencoord=0
alienchrid=1
id1=1
id2=1
fn1=outn+'_native'
fn2=outn+'_alien'
fout1=open(fn1,'w')
fout2=open(fn2,'w')

chrname='scaffold_'

with open(infile) as fin:
	for line in fin:
		lst=line.rstrip().split('\t')
		if len(lst)==1:
			print '{2}{0}:{1}'.format(alienchrid,aliencoord,chrname)
			aliencoord=0
			alienchrid+=1
			continue
		a=int(lst[1])
		b=int(lst[2])

		if a>=b:
			print 'wrong line: '+line
			sys.exit()

		# native
		fout1.write('{0}\t{1}\t{2}\tid:{3},genomealign:{{chr:"{8}{4}",start:{5},stop:{6},strand:"{7}"}}\n'.format(
			lst[0],a,b,
			id1,
			alienchrid,
			aliencoord,
			aliencoord+b-a,
			lst[3],
			chrname
			))
		id1+=1
		# alien
		fout2.write('{8}{0}\t{1}\t{2}\tid:{3},genomealign:{{chr:"{4}",start:{5},stop:{6},strand:"{7}"}}\n'.format(
			alienchrid,
			aliencoord,
			aliencoord+b-a,
			id2,
			lst[0],a,b,
			lst[3],
			chrname
			))
		id2+=1
		aliencoord+=b-a

print '{2}{0}:{1}'.format(alienchrid,aliencoord,chrname)

fout1.close()
fout2.close()

import os

os.system('sort -k1,1 -k2,2n '+fn1+' > x')
os.system('mv x '+fn1)
os.system('bgzip -f '+fn1)
os.system('tabix -f -p bed '+fn1+'.gz')

os.system('sort -k1,1 -k2,2n '+fn2+' > x')
os.system('mv x '+fn2)
os.system('bgzip -f '+fn2)
os.system('tabix -f -p bed '+fn2+'.gz')
