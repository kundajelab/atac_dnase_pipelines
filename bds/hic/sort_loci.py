#!/usr/bin/env python
'''
Created on Aug 08, 2012

@author: ferhat
'''
import sys

USAGE = """USAGE: sort_loci.py <interactionsFile> <numberOfLoci>

Reads the interactions file which contains rows that correspond to interactions between 2 or more loci and sorts these loci according to their chromosome names first and then according to their start indices and then according to their strands. 

<interactionsFile> is the main input file that contains interactions in this format (an example of 3 loci):
"FCD059NACXX:5:2107:18137:185323 +       chr10   100000018       GAA..       +       chr5    66808096	ATG.. -       chr5    808096        ATC.." 

Output will be of the same format but just sorted. For teh above example:
"FCD059NACXX:5:2107:18137:185323 -       chr5    808096        ATC.. +       chr5    66808096    ATG.. +       chr10   100000018       GAA.."

"""


def main(infilename, numberOfLoci):
	if (infilename == "-"):
		infile = sys.stdin
	else:
		infile = open(infilename, "r")
	linecount=0
	for inline in infile:
	#"FCD059NACXX:5:2107:18137:185323 +       chr10   100000018       GAA..       +       chr5    66808096	ATG.. -       chr6    66808096        ATC.."
		if len(inline)<2:
			sys.stderr.write(repr(linecount)+" lines were read from the file. ")
			return
		linecount+=1
		tokenizedStr=inline.rstrip().split()
		readID=tokenizedStr[0]
		strandSigns=[tokenizedStr[i] for i in range(1,4*numberOfLoci,4)]
		tempchrs=[tokenizedStr[i] for i in range(2,4*numberOfLoci,4)]
		readStarts=[int(tokenizedStr[i]) for i in range(3,4*numberOfLoci,4)]
		seqs=[tokenizedStr[i] for i in range(4,4*numberOfLoci+1,4)]

		#print(tempchrs)
		#print(strandSigns)
		#print(readStarts)

		tempchrs,readStarts,strandSigns,seqs=zip(*sorted(zip(tempchrs,readStarts,strandSigns,seqs)))
		#print(tempchrs)
		#print(strandSigns)
		#print(readStarts)
		print(readID+"\t"),
		for i in range(numberOfLoci):
			print(strandSigns[i]+"\t"+tempchrs[i]+"\t"+repr(readStarts[i])+"\t"+seqs[i]),
		print
	sys.stderr.write(repr(linecount)+" lines were read from the file. \n")

if __name__ == "__main__":
	if (len(sys.argv) != 3):
		sys.stderr.write(USAGE)
		sys.exit(1)
	main(str(sys.argv[1]),int(sys.argv[2]))

