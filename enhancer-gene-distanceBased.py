#!/usr/bin/env python
# coding: utf-8

# In[41]:


import pandas as pd
from pybedtools import BedTool
from pybedtools.featurefuncs import TSS
import os


# In[42]:


import sys
args = sys.argv  # a list of the arguments provided (str)
print("processing distanceBased", args)
rawBed= args[1]


# In[43]:


def getAllgeneWithinDistance(enhancer,gene,distance):
    snp = BedTool(enhancer)
    gene = BedTool(gene)
    dfSNP = pd.read_csv(enhancer, sep='\t', header=None)
    win = snp.window(gene,w=10000,output='tempDistance.bed')
    df1 = pd.read_csv('tempDistance.bed', sep='\t', header=None)
    df1.drop(columns=df1.columns[:len(dfSNP.columns)], 
        axis=1, 
        inplace=True)
    df1.to_csv('result.tsv', sep="\t")
    
    os.remove('tempDistance.bed')


# In[44]:


pd.set_option('display.max_columns', None)


# In[45]:


def getClosestgeneWithinDistance(enhancer,genes,distance,num):
    snp = BedTool(enhancer)
    gene = BedTool(genes)
    gene.sort()  
    nearby = snp.closest(gene, d=True,k=num,io=True, output='tempDistance.bed')
    dfSNP = pd.read_csv(enhancer, sep='\t', header=None)
    dfGene = pd.read_csv(genes, sep='\t', header=None)
    df1 = pd.read_csv('tempDistance.bed', sep='\t', header=None)
    df1 = df1[(df1[len(df1.columns)-1].between(-distance, distance))]

    
    df1.drop(columns=df1.columns[4:len(dfSNP.columns)], 
        axis=1, 
        inplace=True)
    df1.drop(columns=df1.columns[8:len(dfGene.columns)+4], 
        axis=1, 
        inplace=True)
    os.remove('tempDistance.bed')
    return df1


# In[46]:


df32 = pd.read_csv('mart_export.txt.gz',sep="\t")


# In[47]:


df32.head()


# In[48]:


def startPosition(df):
    return 'chr'+str(df['Chromosome/scaffold name'])


# In[49]:


def startPosition(df):
    return 'chr'+str(df['Chromosome/scaffold name'])


# In[50]:


df32['0'] = df32.apply(startPosition,axis=1)
df32['1'] = df32['Gene start (bp)']
df32['2'] = df32['Gene end (bp)']
df32['3'] = df32['Gene stable ID']


# In[51]:


df32.tail()


# In[52]:


df21 = df32[['0','1','2','3']]


# In[53]:


df21.head()


# In[54]:


df22 = df21.drop_duplicates(
  subset = ['0','1','2'],
  keep = 'last').reset_index(drop = True)


# In[55]:


df22.tail()


# In[56]:


df21 = df21.sort_values(by=['0','1'],ascending=True)


# In[57]:


new_header = df21.iloc[0] 
df21 = df21[1:] 
df21.columns = new_header 


# In[58]:


df21.to_csv('temp21.bed', sep="\t",index=False)


# In[59]:


#closestGene = getClosestgeneWithinDistance('ENCFF979PSS.bed.gz','../../../storage/data/dna/human/hg38/hg38.trf.bed.gz',1000000,1)


# In[60]:


# closestGene1 = getClosestgeneWithinDistance('ENCFF979PSS.bed.gz','temp21.bed',1000000,1)
closestGene1 = getClosestgeneWithinDistance(rawBed,'temp21.bed',1000000,1)


# In[61]:


# closestGene.to_csv('result.tsv', sep="\t")


# In[62]:


# closestGene.to


# In[63]:


closestGene1.tail()


# In[64]:


df32[(df32['Gene start (bp)'] ==147629652) & (df32['Gene end (bp)'] == 147670524) ]


# ### Duplication issue
# Multiple genes with same gene stable ID , gene start and gene end but different Transcript stable ID

# In[65]:


closestGene1[closestGene1[3].isin(['EH38E2836060','EH38E2836065','EH38E1382733','EH38E2836036','EH38E2836028','EH38E2835994'])]


# In[ ]:





# In[66]:


closestGene1.tail()


# In[67]:


closestGene1[closestGene1[3].isin(['EH38E1310568','EH38E1310568','EH38E1310568','EH38E1310568','EH38E2777292'])]


# In[68]:


closestGene1[closestGene1[3] =='EH38E2789351']


# In[69]:


closestGene1 = closestGene1[closestGene1[14].str.contains('ENSG')]


# In[70]:


closestGene1 = closestGene1[closestGene1[3].str.contains('EH')]


# In[71]:


closestGene1 = closestGene1.drop_duplicates().reset_index(drop = True)


# In[72]:


closestGene1.to_csv('tempFinaldistance.bed',sep='\t',header=None,index=False)


# In[73]:


len(closestGene1.nunique())


# In[74]:


df32 = pd.read_csv('mart_export.txt.gz',sep="\t")


# In[75]:


df32.head()


# In[76]:


df32[(df32['Gene stable ID']== 'ENSG00000162836') ]


# In[77]:


df32[(df32['Gene start (bp)']== 147629652) ]

