#!/usr/bin/env python
# coding: utf-8
import os.path
import uuid
import logging
# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import swifter
import seaborn as sns
import matplotlib.pyplot as plt1
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib_venn import venn2
from collections import Counter

# In[38]:


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
logger = logging.getLogger('waitress')
logger.setLevel(logging.INFO)

# In[39]:
def startPoint(chiapetVal,distanceVal,eqtlVal,imagesFileName):
    logger.info("started visulaization")
    chiaPet = pd.read_csv(chiapetVal, sep='\t', header=None)
    distance = pd.read_csv(distanceVal, sep='\t', header=None)
    eqTL = pd.read_csv(eqtlVal, sep='\t', header=None)
    distance[8] = distance[8].apply(lambda x: np.log10(x) if x != 0 else x)
    chiaPet[9] = 'chiaPet'
    chiaPet[10] = chiaPet.swifter.apply(startPosition, axis=1)
    distance[10] = distance.apply(startPosition, axis=1)
    distance[9] = 'distance'
    eqTL[9] = 'eqtl'
    eqTL[10] = eqTL.apply(startPosition, axis=1)
    finalConcat = pd.concat([chiaPet, eqTL, distance])
    fig1 = sns.countplot(x=9, data=finalConcat)
    fig1.figure.savefig(imagesFileName+"TotalcountComparsion.png")
    fig1.figure.clf()
    finalConcat.groupby([9]).nunique()[3]
    fig27 = finalConcat.groupby([9]).nunique()[3].plot(kind='bar', color="indigo", title='Enhancers',
                                                       ylabel='Unique Enhancer Count',
                                                       xlabel='Method', figsize=(10, 10))
    enhancersFig = fig27.get_figure()
    enhancersFig.suptitle('Enhancer comparision', fontsize=20)

    enhancersFig.savefig(imagesFileName+"enhancers.png")
    fig27.figure.clf()

    fig37 = finalConcat.groupby([9]).nunique()[7].plot(kind='bar', color="red", title='Genes',
                                                       ylabel='Unique Genes Count',
                                                       xlabel='Method', figsize=(10, 10))
    fig37.figure.suptitle('Gene comparision', fontsize=20)

    fig37.figure.savefig(imagesFileName+"GeneComparsion.png")
    fig37.figure.clf()
    ax1 = finalConcat.groupby([9]).nunique()[3].plot(kind='bar', color="indigo", title='Enhancers',
                                                     figsize=(6, 5))
    ax2 = finalConcat.groupby([9]).nunique()[7].plot(kind='bar', color="red", title='Genes',
                                                     figsize=(6, 5))
    indigo = mpatches.Patch(color='indigo', label='Enhancer')

    red = mpatches.Patch(color='red', label='Gene')

    plt.legend(handles=[indigo, red])
    plt.title("Enhancer gene bar")
    # plt.show()
    plt.clf()
    ax1 = finalConcat.groupby([9]).nunique()[3].plot(kind='bar', color="indigo", title='Enhancers', logy=True,
                                                     figsize=(6, 5))
    ax2 = finalConcat.groupby([9]).nunique()[7].plot(kind='bar', color="red", title='Genes', logy=True,
                                                     figsize=(6, 5))
    indigo = mpatches.Patch(color='indigo', label='Enhancer')

    red = mpatches.Patch(color='red', label='Gene')
    plt1.legend(handles=[indigo, red])
    plt1.savefig(imagesFileName+'enhancerGene.png')
    plt1.clf()
    allDF = chiaPet.merge(eqTL, on=3).merge(distance, on=3)
    allDF1 = allDF[(allDF['7_x'] == allDF['7_y']) & (allDF['7_x'] == allDF[7])]
    k = allDF1[['8_x', '8_y']]
    k = k.rename(columns={'8_x': 'chiaPet', '8_y': 'eqtl'})
    set1 = set(chiaPet[7])
    set2 = set(distance[7])
    set3 = set(eqTL[7])
    fig = plt.figure(figsize=(10, 10))

    venn3([set1, set2, set3], ('chiapet', 'distance', 'eqtl'))
    fig.suptitle('Venn diagram of unique Genes', fontsize=20)
    plt.savefig(imagesFileName+'AllGeneComparsion.png')
    plt.clf()

    set1 = set(chiaPet[3])
    set2 = set(distance[3])
    set3 = set(eqTL[3])
    fig = plt.figure(figsize=(10, 10))
    fig.suptitle('Venn diagram of number of unique Enhancer', fontsize=20)

    venn3([set1, set2, set3], ('chiapet', 'distance', 'eqtl'))
    plt.savefig(imagesFileName+'AllEnhancerComparsion.png')
    plt.clf()

    set1 = set(chiaPet[10])
    set3 = set(distance[10])
    set2 = set(eqTL[10])
    sets = [set1, set2, set3]
    fig = plt.figure(figsize=(10, 10))
    fig.suptitle('Venn diagram of number of unique Enhancer-Gene', fontsize=20)
    ax = plt.gca()
    v = venn3([set1, set2, set3], ('chiaPet', 'eqTL', 'distance'), ax=ax)
    h, l = [], []
    for i in ['10', '01', '111']:
        # remove label by setting them to empty string:
        v.get_label_by_id(i).set_text("")
        # append patch to handles list
        h.append(v.get_patch_by_id(i))
        # append count to labels list

    l.append(chiaPet[8].mean())
    l.append(eqTL[8].mean())
    l.append(distance[8].mean())
    # create legend from handles and labels    
    ax.legend(handles=h, labels=l, title="p-value")

    plt.savefig(imagesFileName+'enhancerGeneVenn.png')
    plt.clf()

    set1 = set(chiaPet[10])
    set2 = set(distance[10])
    set3 = set(eqTL[10])
    fig = plt.figure(figsize=(10, 10))
    fig.suptitle('Venn diagram of number of unique Enhancer', fontsize=20)

    ax = venn2([set1, set3], ('chiapet', 'eqtl'))
    # plt.show()
    plt.savefig(imagesFileName+'enhancerchiaPetEQTL.png')
    plt.clf()
    df1 = dict(Counter(chiaPet[3]))
    df2 = dict(Counter(distance[3]))
    df3 = dict(Counter(eqTL[3]))
    newdf = pd.DataFrame.from_dict(df1, orient="index").reset_index()
    newdf = newdf.rename(columns={0: 'repetitions'})
    newdf[1] = 'chiaPet'
    newdf_distance = pd.DataFrame.from_dict(df2, orient="index").reset_index()
    newdf_distance = newdf_distance.rename(columns={0: 'repetitions'})
    newdf_distance[1] = 'distance'
    newdf_eqtl = pd.DataFrame.from_dict(df3, orient="index").reset_index()
    newdf_eqtl = newdf_eqtl.rename(columns={0: 'repetitions'})
    newdf_eqtl[1] = 'eqtl'
    del (df1)
    del (df2)
    del (df3)
    del (chiaPet)
    del (distance)
    del (eqTL)
    fig22 = sns.histplot(data=newdf, x='repetitions')
    chiapetFig = fig22.get_figure()
    chiapetFig.suptitle('Histogram for chiapet of repetitions of Enhancer-Gene', fontsize=20)
    chiapetFig.savefig(imagesFileName+"chiapetHisto.png")
    newdf_distance[newdf_distance.repetitions > 1].count()
    fig223 = sns.histplot(data=newdf_distance, x='repetitions', bins=5)
    distanceFig = fig223.get_figure()
    distanceFig.suptitle('Histogram for Distance of repetitions of Enhancer-Gene', fontsize=20)
    distanceFig.savefig(imagesFileName+"distanceHisto.png")
    fig223.figure.clf()
    fig228 = sns.histplot(data=newdf_eqtl, x='repetitions')
    eqtlFig = fig228.get_figure()
    eqtlFig.suptitle('Histogram for Distance of repetitions of Enhancer-Gene', fontsize=20)

    eqtlFig.savefig(imagesFileName+"eqtlHisto.png")
    fig228.figure.clf()


def startPosition(df):
    return df[3]+'-'+df[7]

def startPosition1(df):
    return df[3]+'-'+df[8]









# ## Questions 
# ### What do we want
# 1. Maximum number of enhancer-gene combination?
# 2. More number of genes to be linked?
# 3. More number of enhancer to be identified?
# 
# ## Questions to think
# 1. We are only dealing with numbers right now.
# 2. What about quality?
# 3. Quantifying only numbers, is it correct? 
# 4. Distance is like a sore thumb here.

# 1. Removed dot from eqtl
# 2. Data change - treatment of negative distance
# 3. Unique enhancer, unique enhancer





# In[83]:





# Enhancer-Gene combination used

# In[84]:


from collections import Counter


# In[85]:





# In[86]:


# sns.set(rc={'figure.figsize':(11.7,8.27)})

# sns.lineplot(data=k)


# In[87]:





# In[88]:






# In[91]:





# In[ ]:





# # P value analysis(common sets)
# # distance as qualitative in distance method
# # P value display on diagram
# # filter more data (remove dots)
