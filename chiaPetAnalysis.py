import logging

import pandas as pd
from pybedtools import BedTool
import os
import uuid
logger = logging.getLogger('waitress')
logger.setLevel(logging.INFO)
def getKey(df):
    return str(df[0])+'_'+str(df[1])+'_'+str(df[2])+'_'+str(df[3])+'_'+str(df[4])+'_'+str(df[5])


def getClosestgeneWithinDistance(enhancer,genes,distance,num):
    snp = BedTool(enhancer)
    gene = BedTool(genes)
    gene.sort()
    tempfileClosest = 'temp/temp'+str(uuid.uuid4())+'.bed'
    nearby = snp.closest(gene, d=True,k=num, output=tempfileClosest)
    dfSNP = pd.read_csv(enhancer, sep='\t', header=None)
    dfGene = pd.read_csv(genes, sep='\t', header=None)
    df1 = pd.read_csv(tempfileClosest, sep='\t', header=None)
    df1 = df1[df1[len(df1.columns)-1] == distance]
    os.remove(tempfileClosest)
    return df1

def checkIfgeneOrEnhancer(file1,file2,final11,final12):
    a_enhancer = getClosestgeneWithinDistance(file1,final11,0,1)
    b_gene = getClosestgeneWithinDistance(file2,final11,0,1)
    c_enhancer = getClosestgeneWithinDistance(file1,final12,0,1)
    d_gene = getClosestgeneWithinDistance(file2,final12,0,1)
    return [a_enhancer,b_gene,c_enhancer,d_gene]



# In[118]:


def startPosition(df):
    return 'chr'+str(df['Chromosome/scaffold name'])


# In[119]:


def getTSS(df):
    return int(df['Gene start (bp)'])-1


# In[120]:
def startPoint(rawBed,file,tempFileName):
    df1 = pd.read_csv('loops-hg38/' + file, sep='\t', header=None)
    df1['key'] = df1.apply(getKey, axis=1)
    df1.index.get_loc(df1.iloc[2].name)
    df2 = df1[[3, 4, 5, 'key', 6]]
    df1 = df1[[0, 1, 2, 'key', 6]]
    df2 = df2.sort_values(by=[3, 4], ascending=True)
    df1 = df1.sort_values(by=[0, 1], ascending=True)
    final11 = 'temp/temp11'+str(uuid.uuid4())+'.bed'
    final12 = 'temp/temp21'+str(uuid.uuid4())+'.bed'
    df2.to_csv(final11, sep='\t', index=False, header=None)
    df1.to_csv(final12, sep='\t', index=False, header=None)
    df32 = pd.read_csv('mart_export.txt.gz', sep='\t')
    df32['0'] = df32.apply(startPosition, axis=1)
    df32['1'] = df32.apply(getTSS, axis=1)
    df32['2'] = df32['Gene start (bp)']
    df32['3'] = df32['Gene end (bp)']
    df32['4'] = df32['Gene stable ID']
    df21 = df32[['0', '1', '2', '3', '4']]
    df21 = df21.sort_values(by=['0', '1'], ascending=True)
    new_header = df21.iloc[0]
    df21 = df21[1:]
    df21.columns = new_header
    df21.to_csv('tempGene.bed', sep="\t", index=False)
    allEnv = checkIfgeneOrEnhancer(rawBed, 'tempGene.bed',final11,final12)
    df1_isEnhancer = allEnv[0]
    df1_isGene = allEnv[1]
    df2_isEnhancer = allEnv[2]
    df2_isGene = allEnv[3]
    df_all = pd.merge(df1_isEnhancer, df2_isGene, how='inner', left_on=[14], right_on=[8])
    df1_all = pd.merge(df2_isEnhancer, df1_isGene, how='inner', left_on=[14], right_on=[8])
    df_final = pd.concat([df_all, df1_all], axis=0)
    k = df_final['9_y']
    df_final.drop(columns=df_final.columns[4:len(df1_isEnhancer.columns)],
                  axis=1,
                  inplace=True)
    df_final.drop(columns=df_final.columns[9:len(df2_isGene.columns) + len(df1_isEnhancer.columns)],
                  axis=1,
                  inplace=True)
    df_final.drop(columns=df_final.columns[5],
                  axis=1,
                  inplace=True)
    df_final[5] = k
    df_final = df_final.drop_duplicates(subset=['3_x', '4_y'], keep='last')
    new_header = df_final.iloc[0]
    df_final = df_final[1:]
    df_final.columns = new_header
    df_final.to_csv(tempFileName, sep='\t', header=None, index=False)
    os.remove(final11)
    os.remove(final12)
    logger.info('finished processing chiaPet with file '+tempFileName)
    return tempFileName








# print(df1_isEnhancer.head())


# In[130]:






# pd.set_option('display.max_columns', None)
# pd.set_option('display.max_rows', None)


# In[132]:





# In[133]:


# len(df1_isEnhancer.columns)
# len(df2_isGene.columns)


# In[134]:





# In[135]:


# df_final.head()


# In[136]:


# df_final = pd.read_csv('tempFinalchiaPet.bed',sep='\t',header=None)
# df_final.tail()
# df_final[df_final[3].isin(['EH38E2836060','EH38E2836065','EH38E1382733','EH38E2836036','EH38E2836028','EH38E2835994'])]

