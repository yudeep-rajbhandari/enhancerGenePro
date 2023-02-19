#!/usr/bin/env python
# coding: utf-8

# In[97]:


import pandas as pd
from pybedtools import BedTool
from pybedtools.featurefuncs import TSS
import os
import numpy as np


# In[98]:


import sys
args = sys.argv  # a list of the arguments provided (str)
print("processing Chiapet", args)
rawBed= args[1]
file = args[2]


# In[99]:


# df1 = pd.read_csv('ENCFF360QPK.bedpe.gz', sep='\t',header=None)
# df1 = pd.read_csv('Leung_2015.Liver.hg38.peakachu-merged.loops', sep='\t',header=None)
df1 = pd.read_csv('loops-hg38/'+file, sep='\t',header=None)


# In[100]:


df1.head()


# In[101]:


def getKey(df):
    return str(df[0])+'_'+str(df[1])+'_'+str(df[2])+'_'+str(df[3])+'_'+str(df[4])+'_'+str(df[5])


# In[102]:


df1.head()


# In[103]:


df1['key'] = df1.apply(getKey,axis=1)


# In[104]:


df1.head()


# In[105]:


df1.index.get_loc(df1.iloc[2].name)


# In[106]:


df2 = df1[[3,4,5,'key',6]]


# In[107]:


df1 = df1[[0,1,2,'key',6]]


# In[108]:


df2 = df2.sort_values(by=[3,4],ascending=True)
df1 = df1.sort_values(by=[0,1],ascending=True)
df2.to_csv('temp11.bed',sep='\t',index=False,header=None)
df1.to_csv('temp21.bed',sep='\t',index=False,header=None)


# In[109]:


df2.head()


# In[110]:


def getClosestgeneWithinDistance(enhancer,genes,distance,num):
    snp = BedTool(enhancer)
    gene = BedTool(genes)
    gene.sort()  
    nearby = snp.closest(gene, d=True,k=num, output='temp.bed')
    dfSNP = pd.read_csv(enhancer, sep='\t', header=None)
    dfGene = pd.read_csv(genes, sep='\t', header=None)
    df1 = pd.read_csv('temp.bed', sep='\t', header=None)
    df1 = df1[df1[len(df1.columns)-1] == distance]
    os.remove('temp.bed')
    return df1


# In[111]:


# aa = pd.read_csv('ENCFF979PSS.bed.gz',sep='\t',header=None)


# In[112]:


# aa.head()


# In[113]:


df32= pd.read_csv('mart_export.txt.gz',sep='\t')


# In[114]:


df32.head()


# In[115]:


def isEnhancer(a,b,c):
    cc= aa[(aa[0] == a) & (aa[1].between(b,c))]
    return cc


# In[116]:


def isGene(a,b,c):
    cc= bb[(bb[0] == a) & bb[1].between(b,c)]
    return cc


# In[117]:


def checkIfgeneOrEnhancer(file1,file2):
    a_enhancer = getClosestgeneWithinDistance(file1,'temp11.bed',0,1)
    b_gene = getClosestgeneWithinDistance(file2,'temp11.bed',0,1)
    c_enhancer = getClosestgeneWithinDistance(file1,'temp21.bed',0,1)
    d_gene = getClosestgeneWithinDistance(file2,'temp21.bed',0,1)
    return [a_enhancer,b_gene,c_enhancer,d_gene]
    


# In[118]:


def startPosition(df):
    return 'chr'+str(df['Chromosome/scaffold name'])


# In[119]:


def getTSS(df):
    return int(df['Gene start (bp)'])-1


# In[120]:


df32['0'] = df32.apply(startPosition,axis=1)
df32['1'] = df32.apply(getTSS,axis=1)
df32['2'] = df32['Gene start (bp)']
df32['3'] = df32['Gene end (bp)']
df32['4'] = df32['Gene stable ID']


# In[121]:


df32.head()


# In[122]:


df21 = df32[['0','1','2','3','4']]


# In[123]:


df21 = df21.sort_values(by=['0','1'],ascending=True)


# In[ ]:





# In[124]:


new_header = df21.iloc[0] 
df21 = df21[1:] 
df21.columns = new_header 


# In[125]:


df21.head()


# In[126]:


df21.to_csv('tempGene.bed', sep="\t",index=False)


# In[ ]:


print('here '+ rawBed)


# In[127]:


allEnv = checkIfgeneOrEnhancer(rawBed,'tempGene.bed')
# allEnv = checkIfgeneOrEnhancer('ENCFF979PSS.bed.gz','tempGene.bed')


# In[ ]:





# In[128]:


df1_isEnhancer = allEnv[0]
df1_isGene = allEnv[1]
df2_isEnhancer = allEnv[2]
df2_isGene = allEnv[3]


# In[129]:


print('here')
print(df1_isEnhancer.head())


# In[130]:


df_all = pd.merge(df1_isEnhancer, df2_isGene,  how='inner', left_on=[14], right_on = [8])
df1_all = pd.merge(df2_isEnhancer, df1_isGene,  how='inner', left_on=[14], right_on = [8])


# In[131]:


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


# In[132]:


df_final = pd.concat([df_all, df1_all], axis=0)


# In[133]:


len(df1_isEnhancer.columns)
len(df2_isGene.columns)


# In[134]:


k = df_final['9_y']


# In[135]:


df_final.head()


# In[136]:


df_final.drop(columns=df_final.columns[4:len(df1_isEnhancer.columns)], 
        axis=1, 
        inplace=True)
df_final.drop(columns=df_final.columns[9:len(df2_isGene.columns)+len(df1_isEnhancer.columns)], 
        axis=1, 
        inplace=True)
df_final.drop(columns=df_final.columns[5], 
        axis=1, 
        inplace=True)


# In[137]:


df_final[5] = k


# In[138]:


df_final = df_final.drop_duplicates(subset=['3_x', '4_y'], keep='last')


# In[146]:





# In[139]:


new_header = df_final.iloc[0] 
df_final = df_final[1:] 
df_final.columns = new_header 


# In[140]:


df_final.head()


# In[141]:


df_final.to_csv('tempFinalchiaPet.bed',sep='\t',header=None,index=False)


# In[142]:


df_final = pd.read_csv('tempFinalchiaPet.bed',sep='\t',header=None)


# In[143]:


df_final.tail()


# In[144]:


df_final[df_final[3].isin(['EH38E2836060','EH38E2836065','EH38E1382733','EH38E2836036','EH38E2836028','EH38E2835994'])]

