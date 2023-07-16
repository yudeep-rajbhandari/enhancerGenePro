#!/usr/bin/env python
# coding: utf-8
import logging
import math
import os

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt1
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2
from matplotlib_venn import venn3

# In[1]:

# In[38]:


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
logger = logging.getLogger('waitress')
logger.setLevel(logging.INFO)

def checkAndSupplyFile(file):
    file1 = pd.DataFrame(columns=[0, 1, 2, 3, 4, 5, 6, 7, 8,9,10])
    try:
        return pd.read_csv(file, sep='\t', header=None)
    except pd.errors.EmptyDataError as e:
        return file1

def createVenn(chiapet,distance,eqtl):
    ax = plt.gca()
    num=3
    setAll = [chiapet,distance,eqtl]
    li = ['Peakachu', 'Distance', 'eQTL']
    if(len(chiapet)==0):
        num=num-1
        setAll.remove(chiapet)
        li.remove('Peakachu')
    if(len(distance)==0):
        num=num-1
        setAll.remove(distance)
        li.remove('Distance')
    if(len(eqtl)==0):
        num=num-1
        setAll.remove(eqtl)
        li.remove('eQTL')
    if(num==3):
        return venn3(setAll, tuple(li),ax=ax)
    elif (num==2):
        return venn2(setAll,tuple(li),ax=ax)
# In[39]:
def startPoint(chiapetVal,distanceVal,eqtlVal,imagesFileName):
    logger.info("started visulaization")
    chiaPet = checkAndSupplyFile(chiapetVal)
    distance = checkAndSupplyFile(distanceVal)
    eqTL = checkAndSupplyFile(eqtlVal)
    distance[8] = distance[8].apply(lambda x: np.log10(x) if x != 0 else x)
    chiaPet[9] = 'Peakachu'
    if (len(chiaPet[3]) > 0) and (len(chiaPet[7]) > 0):
        chiaPet[10] = chiaPet.apply(startPosition, axis=1)
    if (len(distance[3]) > 0) and (len(distance[7]) > 0):
        distance[10] = distance.apply(startPosition, axis=1)
    distance[9] = 'Distance'
    eqTL[9] = 'eQTL'
    imagesFileNameBar = imagesFileName+'BarPlots/'
    if not os.path.exists(imagesFileNameBar):
        # If it doesn't exist, create it
        os.makedirs(imagesFileNameBar)
    if (len(eqTL[3]) > 0) and (len(eqTL[7]) > 0):
        eqTL[10] = eqTL.apply(startPosition, axis=1)
    finalConcat = pd.concat([chiaPet, eqTL, distance])
    fig1 = sns.countplot(x=9, data=finalConcat)
    fig1.set(title='Enhancer-Gene Links Counts Comparison',xlabel='',ylabel='Count')
    fig1.figure.savefig(imagesFileNameBar+"TotalcountComparsion.pdf",format="pdf", bbox_inches="tight")
    fig1.figure.clf()
    finalConcat.groupby([9]).nunique()[3]
    fig27 = finalConcat.groupby([9]).nunique()[3].plot(kind='bar', color="indigo", title='Enhancers',
                                                       ylabel='# Unique Enhancers',
                                                       xlabel='', figsize=(10, 10))
    enhancersFig = fig27.get_figure()
    # enhancersFig.suptitle('Enhancer comparision', fontsize=20)
    enhancersFig.figure.suptitle('Enhancer comparision', fontsize=20)
    enhancersFig.savefig(imagesFileNameBar+"uniqueEnhancersByMethod.pdf",format="pdf", bbox_inches="tight")
    fig27.figure.clf()

    fig37 = finalConcat.groupby([9]).nunique()[7].plot(kind='bar', color="red", title='Genes',
                                                       ylabel='# Unique Genes',
                                                       xlabel='', figsize=(10, 10))
    fig37.figure.suptitle('Gene comparision', fontsize=20)

    fig37.figure.savefig(imagesFileNameBar+"uniqueGenesByMethod.pdf",format="pdf", bbox_inches="tight")
    fig37.figure.clf()
    ax1 = finalConcat.groupby([9]).nunique()[3].plot(kind='bar', color="indigo", title='Enhancers',
                                                     figsize=(6, 5))
    ax2 = finalConcat.groupby([9]).nunique()[7].plot(kind='bar', color="red", title='Genes',
                                                     figsize=(6, 5))
    indigo = mpatches.Patch(color='indigo', label='Enhancer')

    red = mpatches.Patch(color='red', label='Gene')

    plt.legend(handles=[indigo, red])
    # plt.suptitle("Enhancer gene bar")
    # plt.show()
    plt.clf()
    ax1 = finalConcat.groupby([9]).nunique()[3].plot(kind='bar', color="indigo", title='Enhancers', logy=True,
                                                     figsize=(6, 5))
    ax2 = finalConcat.groupby([9]).nunique()[7].plot(kind='bar', color="red", title='Unique Enhancers Vs Genes', logy=True,
                                                     figsize=(6, 5))
    indigo = mpatches.Patch(color='indigo', label='Enhancer')

    red = mpatches.Patch(color='red', label='Gene')
    plt1.legend(handles=[indigo, red])
    plt1.ylabel('')
    plt1.xlabel('')
    # plt1.suptitle("Enhancer gene bars")
    plt1.savefig(imagesFileNameBar+'enhancerGene.pdf',format="pdf", bbox_inches="tight")
    plt1.clf()
    # allDF = chiaPet.merge(eqTL, on=3).merge(distance, on=3)
    # allDF1 = allDF[(allDF['7_x'] == allDF['7_y']) & (allDF['7_x'] == allDF[7])]
    # k = allDF1[['8_x', '8_y']]
    # k = k.rename(columns={'8_x': 'chiaPet', '8_y': 'eqtl'})
    imagesFileNameVenn = imagesFileName+'VennDiagram/'
    if not os.path.exists(imagesFileNameVenn):
        # If it doesn't exist, create it
        os.makedirs(imagesFileNameVenn)
    set1 = set(chiaPet[7])
    set2 = set(distance[7])
    set3 = set(eqTL[7])
    fig = plt.figure(figsize=(10, 10))

    createVenn(set1, set2, set3)
    fig.suptitle('Overlaps between unique genes linked per method', fontsize=20)
    plt.savefig(imagesFileNameVenn+'AllGeneComparsion.pdf',format="pdf", bbox_inches="tight")
    plt.clf()

    set1 = set(chiaPet[3])
    set2 = set(distance[3])
    set3 = set(eqTL[3])
    fig = plt.figure(figsize=(10, 10))
    fig.suptitle('Overlaps between unique enhancers linked per method', fontsize=20)

    createVenn(set1,set2,set3)
    plt.savefig(imagesFileNameVenn+'AllEnhancerComparsion.pdf',format="pdf", bbox_inches="tight")
    plt.clf()

    set1 = set(chiaPet[10])
    set3 =  set(eqTL[10])
    set2 = set(distance[10])
    sets = [set1, set2, set3]
    fig = plt.figure(figsize=(10, 10))
    fig.suptitle('Overlaps between unique enhancers-Gene linkage per method', fontsize=20)
    ax = plt.gca()
    v = createVenn(set1,set2,set3)
    h, l = [], []
    if len(set1) != 0 and len(set2) !=0 and len(set3) !=0:
        labels = ['100','010', '001']
    else:
        labels = ['100', '010']

    for i in labels:
        # remove label by setting them to empty string:
        # if v.get_label_by_id(i) is not None:
        #     v.get_label_by_id(i).set_text(i)
        # append patch to handles list
        h.append(v.get_patch_by_id(i))
        # append count to labels list
    if not math.isnan(chiaPet[8].mean()):
        l.append(0 if math.isnan(chiaPet[8].mean()) else str("{:.3f}".format(chiaPet[8].mean()))+' cScore')
    if not math.isnan(distance[8].mean()):
        l.append(0 if math.isnan(distance[8].mean()) else str("{:.3f}".format(distance[8].mean()))+' log(distance)')
    if not math.isnan(eqTL[8].mean()):
        l.append(0 if math.isnan(eqTL[8].mean()) else str('{:0.3e}'.format(eqTL[8].mean()))+' pVal(eQTL)')
    # create legend from handles and labels    
    ax.legend(handles=h, labels=l, title="Average measures of strength",loc='upper center', bbox_to_anchor=(0.5, -0.05))
    plt.rcParams.update({'font.size': 22})

    plt.savefig(imagesFileNameVenn+'enhancerGeneVenn.pdf',format="pdf", bbox_inches="tight")
    plt.clf()

    # set1 = set(chiaPet[10])
    # set2 = set(distance[10])
    # set3 = set(eqTL[10])
    # fig = plt.figure(figsize=(10, 10))
    # fig.suptitle('Venn diagram of number of unique Enhancer', fontsize=20)
    #
    # ax = venn2([set1, set3], ('Peakachu', 'eQTL'))
    # # plt.show()
    # plt.savefig(imagesFileNameVenn+'enhancerchiaPetEQTL.pdf',format="pdf", bbox_inches="tight")
    # plt.clf()
    df1 = dict(Counter(chiaPet[3]))
    df2 = dict(Counter(distance[3]))
    df3 = dict(Counter(eqTL[3]))
    newdf = pd.DataFrame.from_dict(df1, orient="index").reset_index()
    newdf = newdf.rename(columns={0: 'Number of enhancers per gene'})
    newdf[1] = 'chiaPet'
    newdf_distance = pd.DataFrame.from_dict(df2, orient="index").reset_index()
    newdf_distance = newdf_distance.rename(columns={0: 'Number of enhancers per gene'})
    newdf_distance[1] = 'distance'
    newdf_eqtl = pd.DataFrame.from_dict(df3, orient="index").reset_index()
    newdf_eqtl = newdf_eqtl.rename(columns={0: 'Number of enhancers per gene'})
    newdf_eqtl[1] = 'eqtl'
    del (df1)
    del (df2)
    del (df3)
    del (chiaPet)
    del (distance)
    del (eqTL)
    imagesFileNameHisto = imagesFileName+'histograms/'
    if not os.path.exists(imagesFileNameHisto):
        # If it doesn't exist, create it
        os.makedirs(imagesFileNameHisto)
    if newdf.shape[0] > 0:
        fig22 = sns.histplot(data=newdf, x='Number of enhancers per gene')
        chiapetFig = fig22.get_figure()
        chiapetFig.suptitle('Distribution of enhancer number per gene (Peakachu)', fontsize=20)
        chiapetFig.savefig(imagesFileNameHisto+"PeakachuHisto.pdf",format="pdf", bbox_inches="tight")
        chiapetFig.figure.clf()
    # newdf_distance[newdf_distance.repetitions > 1].count()
    if newdf_distance.shape[0] > 0:
        fig223 = sns.histplot(data=newdf_distance, x='Number of enhancers per gene')
        distanceFig = fig223.get_figure()
        distanceFig.suptitle('Distribution of enhancer number per gene (Distance)', fontsize=20)
        distanceFig.savefig(imagesFileNameHisto+"DistanceHisto.pdf",format="pdf", bbox_inches="tight")
        fig223.figure.clf()
    if newdf_eqtl.shape[0] > 0:
        fig228 = sns.histplot(data=newdf_eqtl, x='Number of enhancers per gene')
        eqtlFig = fig228.get_figure()
        eqtlFig.suptitle('Distribution of enhancer number per gene (eQTL)', fontsize=20)
        eqtlFig.savefig(imagesFileNameHisto+"eqtlHisto.pdf",format="pdf", bbox_inches="tight")
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
