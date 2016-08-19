#!/usr/bin/env python2

import sys,os
sys.path.append('/home/xzhou/subtleKnife/script/genescript')
import parseUcscgenestruct

if len(sys.argv)!=3:
	print '<ensGene.txt file> <tkname> knownToEnsembl.txt and kgXref.txt must be under current dir'
	sys.exit()


aa={}
with open('knownToEnsembl.txt') as fin:
	for line in fin:
		lst=line.rstrip().split('\t')
		aa[lst[0]]=lst[1]

symbol={}
desc={}
with open('kgXref.txt') as fin:
	for line in fin:
		lst=line.rstrip().split('\t')
		if lst[0] in aa:
			ens=aa[lst[0]]
			if len(lst[4])>0:
				symbol[ens]=lst[4]
			if len(lst[7])>0:
				desc[ens]=lst[7]


ucsc,tkname=sys.argv[1:]



# dump
fout=open(tkname,'w')
fout2=open(tkname+'_load','w')

id=1
with open(ucsc) as fin:
	for line in fin:
		lst=line.rstrip().split('\t')
		g=parseUcscgenestruct.parse(lst,True)
		name=lst[1]
		fout.write('{0}\t{1}\t{2}\tname:"{3}",id:{4},strand:"{5}",'.format(
			g['chrom'],
			g['start'],
			g['stop'],
			name,
			id,
			g['strand']))
		id+=1
		if 'thin' in g or 'thick' in g:
			fout.write('struct:{')
			if 'thin' in g:
				fout.write('thin:[')
				for x in g['thin']:
					fout.write('[{0},{1}],'.format(x[0],x[1]))
				fout.write('],')
			if 'thick' in g:
				fout.write('thick:[')
				for x in g['thick']:
					fout.write('[{0},{1}],'.format(x[0],x[1]))
				fout.write('],')
			fout.write('},')
		# desc
		if name in desc:
			fout.write('desc:"'+desc[name]+'",')
		if name in symbol:
			fout.write('name2:"'+symbol[name]+'"')
			fout2.write('{0}\t{1}\t{2}\t{3}\n'.format(g['chrom'],g['start'],g['stop'],symbol[name]))
		fout.write('\n')
		fout2.write('{0}\t{1}\t{2}\t{3}\n'.format(g['chrom'],g['start'],g['stop'],name))


fout2.close()
fout.close()

import os
os.system('sort -k1,1 -k2,2n '+tkname+' > x')
os.system('mv x '+tkname)
os.system('bgzip -f '+tkname)
os.system('tabix -f -p bed '+tkname+'.gz')

print '''
drop table if exists {0};
create table {0} (
chrom varchar(20) not null,
start int unsigned not null,
stop int unsigned not null,
name varchar(100) not null
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
load data local infile '{0}_load' into table {0};
create index name on {0} (name);
'''.format(tkname)

