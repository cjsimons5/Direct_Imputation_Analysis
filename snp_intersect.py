#!/usr/bin/env python

import pandas as pd
import sys

if len(sys.argv)<=4:
    print('python snp_intersect.py file1 file2 chromosome ofile_prefix')
    exit(1)

ImpFile=sys.argv[1]
SeqFile=sys.argv[2]
chrom=sys.argv[3]

#output file prefix should have sequencing data name and imputation data shorthand names to identify which files are being intersected
#    i.e. SIGMA OMNI imputed data and SIGMA 52k sequence data can have prefix OMNI_52k
prefix=sys.argv[4]

Imp=open(ImpFile)
Imp=Imp.read()
ImpSNPList=Imp.splitlines()

Seq=open(SeqFile)
Seq=Seq.read()
SeqSNPList=Seq.splitlines()

def intersect_snps(L1, L2):
    temp = set(L2)
    L3 = [value for value in L1 if value in temp]
    return L3

ofile=open('{}_Keep_chr{}.txt'.format(prefix, chrom), 'w+')
keepers=intersect_snps(ImpSNPList, SeqSNPList)
for snp in keepers:
    ofile.write(snp + '\n')
ofile.close()

logfile=open('{}_Keep_chr{}.log'.format(prefix, chrom), 'w+')
logfile.write('Intersection: ')
logfile.write(str(len(keepers)))
logfile.write('\n')
logfile.write('Imputed SNP total: ')
logfile.write(str(len(ImpSNPList)))
logfile.write('\n')
logfile.write('Sequence SNP total: ')
logfile.write(str(len(SeqSNPList)))
logfile.write('\n')
logfile.close()
