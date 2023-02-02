#!/usr/bin/env python

import pandas as pd
import sys

if len(sys.argv)<=4:
    print('python intersect_names.py file1 file2 outputname1 outputname2')
    exit(1)
    
file1=sys.argv[1]
file2=sys.argv[2]
name1=sys.argv[3]
name2=sys.argv[4]
    
def intersect_names(L1, L2):
    temp = set(L2)
    L3 = [value for value in L1 if value in temp]
    return L3

list1=pd.read_csv(file1, sep='\t')
list2=pd.read_csv(file2, sep='\t')

name_keepers=intersect_names(list2, list1)
ofile=open('intersected_names_{}_{}.txt'.format(name1, name2), 'w')
for name in name_keepers[9:]:
    ofile.write(name+ '\n')
ofile.close()
