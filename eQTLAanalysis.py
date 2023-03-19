#!/usr/bin/env python
# coding: utf-8
import logging
import uuid

# In[1]:


import pandas as pd

from pybedtools import BedTool
from pybedtools.featurefuncs import TSS
import os


# In[1]:


import sys



# In[2]:
def startPoint(rawBed,file,fileHelp,finalTempEQTL):
    df1 = pd.read_csv('GTEx_Analysis_v8_eQTL/'+file, sep='\t')
    df1['start'] = df1.apply(startPosition, axis=1)
    df1['end'] = df1.apply(endPosition, axis=1)
    df1['chrom'] = df1.apply(chrom, axis=1)
    df2 = df1[['chrom', 'start', 'end', 'gene_id', 'pval_nominal']]
    df2 = df2.sort_values(by=['chrom', 'end'], ascending=True)
    new_header = df2.iloc[0]
    df2 = df2[1:]
    df2.columns = new_header
    tempEQTL1 = 'temp1'+str(uuid.uuid4())+'.bed'
    df2.to_csv(tempEQTL1, sep="\t", index=False)
    allBed = getClosestgeneWithinDistance(tempEQTL1, rawBed, 1000000, 1)
    allBed = allBed[allBed[16] <= 100]
    allBed = allBed.sort_values(by=[15, 0], ascending=True)
    df4 = pd.read_csv('GTEx_Analysis_v8_eQTL/' + fileHelp, sep='\t')
    finaldf = allBed.merge(df4, left_on=3, right_on='gene_id', how='left')
    finaldf = finaldf[[4, 5, 6, 7, 8, 'gene_chr', 'gene_start', 'gene_end', 'gene_id', 15]]
    finaldf[15] = finaldf[4]
    finaldf.head()
    finaldf = finaldf.rename(
        columns={5: "enahancer_chr", 6: "enhancer_start", 7: "enhancer_end", 8: "enhancer_id", 15: "p-val"})
    finaldf.drop(columns=finaldf.columns[0],
                 axis=1,
                 inplace=True)
    finaldf['gene_id'] = finaldf.apply(geneIDNorm, axis=1)
    finaldf = finaldf.drop_duplicates(subset=['enhancer_id', 'gene_id'], keep='last')
    finaldf.to_csv(finalTempEQTL, sep='\t', header=None, index=False)
    os.remove(tempEQTL1)
    logging.info("Finished processing eqtl")
    return finalTempEQTL


def startPosition(df):
    return df.variant_id.split('_')[1]


# In[5]:


def endPosition(df):
    return int(df.start)+1


# In[6]:


def chrom(df):
    return df.variant_id.split('_')[0]


# In[7]:


def geneIDNorm(df):
    return df.gene_id.split('.')[0]




def getClosestgeneWithinDistance(enhancer,genes,distance,num):
    snp = BedTool(enhancer)
    gene = BedTool(genes)
    gene.sort()
    tempEQTL = 'tempEQTL'+str(uuid.uuid4())+'.bed'
    nearby = snp.closest(gene, d=True,k=num, output=tempEQTL)
    dfSNP = pd.read_csv(enhancer, sep='\t', header=None)
    dfGene = pd.read_csv(genes, sep='\t', header=None)
    df1 = pd.read_csv(tempEQTL, sep='\t', header=None)
    df1 = df1[(df1[len(df1.columns)-1].between(-distance, distance))]
    
    df1.drop(columns=df1.columns[5:len(dfSNP.columns)], 
        axis=1, 
        inplace=True)
    print(df1.head())
    df1.drop(columns=df1.columns[9:len(dfGene.columns)+4], 
        axis=1, 
        inplace=True)
    os.remove(tempEQTL)
    return df1


# In[18]:



