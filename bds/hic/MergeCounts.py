#!/usr/bin/env python

# Combines FixedWindowsMidpoints and ContactsPerFragment. For windows that are not in ContactsPerFragment, interaction count is assigned to 0. 
import sys
from operator import itemgetter
import csv

USAGE = """USAGE: """

def main(infilename1, infilename2, outfile):
    tokenizedStr1=[]
    infile1 = open(infilename1, "r")
    for inline1 in infile1:
        tokenizedStr1.append(inline1.rstrip().split())

    tokenizedStr2=[]
    infile2 = open(infilename2, "r")
    for inline2 in infile2:
        tokenizedStr2.append(inline2.rstrip().split())
   
    temp_list=[]  
    for i in tokenizedStr2:
        temp_list.append([i[0], i[1]])
#    print temp_list

    new_list=[]
    for i in tokenizedStr1:
        if i not in temp_list:
            new_list.append(i+['0'])
#    print new_list
    
    sorted_list = sorted(tokenizedStr2+new_list, key=itemgetter(0,1))
#    print sorted_list 

    with open(outfile, 'a') as outcsv:   
#configure writer to write standart csv file
        writer = csv.writer(outcsv, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
        for item in sorted_list:
#Write item to outcsv
#        print -(len(item[0]-3)
#            writer.writerow([item[0][-(len(item[0])-3):], '1',  item[1],  item[2], '1']) 
            writer.writerow([item[0], '1',  item[1],  item[2], '1']) 
     
 


if __name__ == "__main__":
        if (len(sys.argv) != 4):
                sys.stderr.write(USAGE)
                sys.exit(1)
        main(str(sys.argv[1]),str(sys.argv[2]), str(sys.argv[3]))
