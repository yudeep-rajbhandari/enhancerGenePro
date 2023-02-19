#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from pybedtools import BedTool
from pybedtools.featurefuncs import TSS
import os


# In[1]:


import sys
args = sys.argv  # a list of the arguments provided (str)
print("processing eqtl", args)
rawBed= args[1]
file = args[2]
fileHelp = args[3]


# In[2]:


df1 = pd.read_csv('GTEx_Analysis_v8_eQTL/'+file, sep='\t')


# In[3]:


df1.head()


# In[4]:


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


# In[8]:


df1['start'] = df1.apply(startPosition,axis=1)
df1['end'] = df1.apply(endPosition,axis = 1)
df1['chrom'] = df1.apply(chrom,axis=1)


# In[ ]:





# In[9]:


df2 = df1[['chrom','start','end','gene_id','pval_nominal']]


# In[ ]:





# In[10]:


df2.head()


# In[11]:


df2 = df2.sort_values(by=['chrom','end'],ascending=True)


# In[12]:


new_header = df2.iloc[0] 
df2 = df2[1:] 
df2.columns = new_header 


# In[13]:


df2.head(10)


# In[14]:


df2.to_csv('temp1.bed', sep="\t",index=False)


# In[15]:


df3 = pd.read_csv('temp1.bed',sep = '\t')


# In[16]:


df3.head()


# In[17]:


def getClosestgeneWithinDistance(enhancer,genes,distance,num):
    snp = BedTool(enhancer)
    gene = BedTool(genes)
    gene.sort()  
    nearby = snp.closest(gene, d=True,k=num, output='tempEQTL.bed')
    dfSNP = pd.read_csv(enhancer, sep='\t', header=None)
    dfGene = pd.read_csv(genes, sep='\t', header=None)
    df1 = pd.read_csv('tempEQTL.bed', sep='\t', header=None)
    df1 = df1[(df1[len(df1.columns)-1].between(-distance, distance))]
    
    df1.drop(columns=df1.columns[5:len(dfSNP.columns)], 
        axis=1, 
        inplace=True)
    print(df1.head())
    df1.drop(columns=df1.columns[9:len(dfGene.columns)+4], 
        axis=1, 
        inplace=True)
    os.remove('tempEQTL.bed')
    return df1


# In[18]:


allBed = getClosestgeneWithinDistance('temp1.bed',rawBed,1000000,1)


# In[19]:


allBed.head()


# In[20]:


len(allBed)


# In[ ]:





# In[21]:


allBed = allBed[allBed[16] <= 100 ]


# In[22]:


allBed.head()


# In[23]:


len(allBed)


# In[24]:


allBed = allBed.sort_values(by=[15,0],ascending=True)


# In[25]:


allBed.head()


# In[26]:


df4 = pd.read_csv('GTEx_Analysis_v8_eQTL/'+fileHelp,sep = '\t')


# In[27]:


df4[(df4['gene_chr'] =='chr1') & (df4['gene_start'] ==147632194) & (df4['gene_end']==147632224)]


# In[28]:


df4.head()


# In[29]:


finaldf = allBed.merge(df4,left_on=3,right_on='gene_id',how='left')


# In[30]:


finaldf.head()


# In[31]:


finaldf = finaldf[[4,5,6,7,8,'gene_chr','gene_start','gene_end','gene_id',15]]


# In[32]:


finaldf.head(10)


# In[33]:


finaldf[15] = finaldf[4]


# In[34]:


finaldf.head()


# In[35]:


finaldf = finaldf.rename(columns={5: "enahancer_chr", 6: "enhancer_start",7:"enhancer_end",8:"enhancer_id",15:"p-val"})


# In[36]:


finaldf.drop(columns=finaldf.columns[0], 
        axis=1, 
        inplace=True)


# In[37]:


finaldf.head()


# ### Duplications issue
# Here the duplication is because more than one eQtl exists in between enhancer start and end  pointing to the same gene in most occasions.
# 
# Solution: enahncerID and geneID should be unique.

# In[38]:


finaldf['gene_id'] = finaldf.apply(geneIDNorm,axis=1)


# In[39]:


finaldf = finaldf.drop_duplicates(subset=['enhancer_id', 'gene_id'], keep='last')


# In[40]:


finaldf[finaldf['enhancer_id'].isin(['EH38E2836060','EH38E2836065','EH38E1382733','EH38E2836036','EH38E2836028','EH38E2835994'])]


# In[41]:


finaldf.tail()


# In[42]:


finaldf.to_csv('tempFinaleQTL.bed',sep='\t',header=None,index=False)

