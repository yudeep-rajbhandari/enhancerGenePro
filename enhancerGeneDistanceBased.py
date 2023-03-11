#!/usr/bin/env python
# coding: utf-8

# In[41]:


import pandas as pd
from pybedtools import BedTool
from pybedtools.featurefuncs import TSS
import os


# In[42]:


import sys
import uuid


def getClosestgeneWithinDistance(enhancer,genes,distance,num):
    snp = BedTool(enhancer)
    gene = BedTool(genes)
    gene.sort()
    tempFileName = 'tempDistance'+str(uuid.uuid4())+'.bed'
    nearby = snp.closest(gene, d=True,k=num,io=True, output=tempFileName)
    dfSNP = pd.read_csv(enhancer, sep='\t', header=None)
    dfGene = pd.read_csv(genes, sep='\t', header=None)
    df1 = pd.read_csv(tempFileName, sep='\t', header=None)
    df1 = df1[(df1[len(df1.columns)-1].between(-distance, distance))]

    
    df1.drop(columns=df1.columns[4:len(dfSNP.columns)], 
        axis=1, 
        inplace=True)
    df1.drop(columns=df1.columns[8:len(dfGene.columns)+4], 
        axis=1, 
        inplace=True)
    os.remove(tempFileName)
    return df1


# In[46]:



def startPosition(df):
    return 'chr'+str(df['Chromosome/scaffold name'])


# In[49]:


def startPosition(df):
    return 'chr'+str(df['Chromosome/scaffold name'])

def startPoint(rawBed,tempFinalFile):
    # df32 = pd.read_csv('mart_export.txt.gz', sep="\t")
    # df32['0'] = df32.apply(startPosition, axis=1)
    # df32['1'] = df32['Gene start (bp)']
    # df32['2'] = df32['Gene end (bp)']
    # df32['3'] = df32['Gene stable ID']
    # df21 = df32[['0', '1', '2', '3']]
    # df22 = df21.drop_duplicates(
    #     subset=['0', '1', '2'],
    #     keep='last').reset_index(drop=True)
    # df21 = df21.sort_values(by=['0', '1'], ascending=True)
    # new_header = df21.iloc[0]
    # df21 = df21[1:]
    # df21.columns = new_header
    # df21.to_csv('tempDistanceBased.bed', sep="\t", index=False)
    closestGene1 = getClosestgeneWithinDistance(rawBed, 'tempDistanceBased.bed', 1000000, 1)
    closestGene1 = closestGene1[closestGene1[14].str.contains('ENSG')]
    closestGene1 = closestGene1[closestGene1[3].str.contains('EH')]
    closestGene1 = closestGene1.drop_duplicates().reset_index(drop=True)
    closestGene1.to_csv(tempFinalFile, sep='\t', header=None, index=False)
    print('finished processing distance based')
    return tempFinalFile






