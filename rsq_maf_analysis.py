#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys
import os

def issame(df):
    for i in range(0, len(df)):
        if df.iloc[i]!=df.iloc[0]:
            return False
    return True

##rsquared calculation##
def getrsq(seq, imp):
    #order of true genotypes file first, imputation file second is important
    rsq=[]
    for i in range(0, seq.shape[0]):
        #remove missing genotype individuals from true genotypes for comparison at SNP i
        miss = seq.columns[seq.iloc[i,:]=='NaN']
        SNP_i_imp = imp.iloc[i,:]
        #remove for imputation file as well, can't compare what does not exist
        SNP_i_imp = SNP_i_imp.drop(miss)
        SNP_i_seq = seq.iloc[i,:]
        SNP_i_seq = SNP_i_seq.drop(miss)
        SNP_i_imp  = pd.to_numeric(SNP_i_imp)
        SNP_i_seq = pd.to_numeric(SNP_i_seq)
        
        #if snp is monomorphic, change the value by 1e-6 to provide sample variance
        if issame(SNP_i_imp):
            if SNP_i_imp[0]==0.0:
                SNP_i_imp[0]=0.000001
            elif SNP_i_imp[0]==1:
                SNP_i_imp[0]=1.000001
            elif SNP_i_imp[0]==2:
                SNP_i_imp[0]=1.999999
                
        #pearson r, and then r squared appended to output list
        corr = np.corrcoef(SNP_i_seq, SNP_i_imp)
        rsq.append(np.power(corr[0][1], 2))
    return(rsq)
##MAF calculation##
def MAF(seq):
    maf=[]
    for i in range(0, seq.shape[0]):
        miss = seq.columns[seq.iloc[i,:]=='NaN']
        SNP_i_seq = seq.iloc[i,:]
        SNP_i_seq = SNP_i_seq.drop(miss)
        app=sum(SNP_i_seq)
        maf.append(app/(2*len(SNP_i_seq)))
    return(maf)

if len(sys.argv)<=4:
    print('python rsq_maf_analysis.py IMPUTED.vcf SEQUENCES.vcf ImpName SeqName')
    exit(1)
ImpFile=sys.argv[1]
SeqFile=sys.argv[2]

#IMPUTED DOSAGE DF
imp = pd.read_csv(ImpFile, sep='\t')
imp = imp.set_index(imp.iloc[:,2])
imp = imp.drop(imp.iloc[:,0:9], axis=1)
imp = imp.reindex(sorted(imp.columns), axis=1)
for column in imp:
    imp[column] = imp[column].str.split(':').str[1]
    imp[column] = pd.to_numeric(imp[column])
    
#SEQUENCES DF
seq = pd.read_csv(SeqFile, sep='\t')
seq = seq.set_index(seq.iloc[:,2])
seq = seq.drop(seq.iloc[:,0:9], axis=1)
seq = seq.reindex(sorted(seq.columns), axis=1)
seq=seq.replace('0/0', 0, regex=True)
seq=seq.replace('1/0', 1, regex=True)
seq=seq.replace('0/1', 1, regex=True)
seq=seq.replace('1/1', 2, regex=True)
seq=seq.replace('./.', "NaN", regex=True)
seq=seq.replace('0\|0', 0, regex=True)
seq=seq.replace('1\|0', 1, regex=True)
seq=seq.replace('0\|1', 1, regex=True)
seq=seq.replace('1\|1', 2, regex=True)
seq=seq.replace('.\|.', "NaN", regex=True)

#Run Analysis
maf=MAF(seq)
rsq=getrsq(seq, imp)
ofile=open('maf_rsq_{}_{}.txt'.format(sys.argv[3], sys.argv[4]), 'a')
fline=open('maf_rsq_{}_{}.txt'.format(sys.argv[3], sys.argv[4]), 'r')
fline=fline.readline()
title='SNP' + '\t' + 'MAF' + '\t' + 'r^2' + '\n'
if fline != title:
    ofile.write('SNP')
    ofile.write('\t')
    ofile.write('MAF')
    ofile.write('\t')
    ofile.write('r^2')
    ofile.write('\n')
for i in range(0, seq.shape[0]):
    ofile.write(seq.index[i])
    ofile.write('\t')
    ofile.write(str(round(maf[i], 5)))
    ofile.write('\t')
    ofile.write(str(round(rsq[i], 5)))
    ofile.write('\n')
ofile.close()
