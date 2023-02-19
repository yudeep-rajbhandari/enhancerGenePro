#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import swifter


# In[38]:


pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)


# In[39]:


chiaPet = pd.read_csv('tempFinalchiaPet.bed', sep='\t',header=None)


# In[40]:


distance = pd.read_csv('tempFinaldistance.bed', sep='\t',header=None)


# In[41]:


eqTL = pd.read_csv('tempFinaleQTL.bed', sep='\t',header=None)


# In[42]:


distance.head(10)


# In[43]:


distance[8] = distance[8].apply(lambda x: np.log10(x) if x != 0 else x)


# In[44]:


def startPosition(df):
    return df[3]+'-'+df[7]


# In[45]:


def startPosition1(df):
    return df[3]+'-'+df[8]


# In[46]:


eqTL.head()


# In[47]:


chiaPet.head()


# In[48]:


# chiaPet.drop(columns=6,
#         axis=1, 
#         inplace=True)


# In[49]:


chiaPet.head()


# In[50]:


chiaPet[9]= 'chiaPet'


# In[51]:


chiaPet[10] = chiaPet.swifter.apply(startPosition,axis=1)


# In[52]:


chiaPet.head()


# In[53]:


len(chiaPet[7].unique())


# In[54]:


len(chiaPet[3].unique())


# In[ ]:





# In[55]:


distance[10] = distance.apply(startPosition,axis=1)


# In[56]:


distance[9] = 'distance'


# In[57]:


len(distance[3].unique())


# In[58]:


len(distance[6].unique())


# In[59]:


eqTL[9] = 'eqtl'


# In[60]:


eqTL[10] = eqTL.apply(startPosition,axis=1)


# In[61]:


eqTL.head()


# In[62]:


import seaborn as sns


# In[63]:


chiaPet.head()


# In[64]:


finalConcat = pd.concat([chiaPet,eqTL,distance])


# In[65]:


finalConcat.head()


# In[66]:


import seaborn as sns


# In[67]:


fig1 = sns.countplot(x =9, data = finalConcat)
fig1.figure.savefig("images/TotalcountComparsion.png") 
fig1.figure.clf()


# In[68]:


finalConcat.groupby([9]).nunique()[3]


# In[69]:


fig27 = finalConcat.groupby([9]).nunique()[3].plot(kind='bar' , color="indigo",title='Enhancers', ylabel='Unique Enhancer Count',
         xlabel='Method', figsize=(10, 10))
enhancersFig = fig27.get_figure()
enhancersFig.suptitle('Enhancer comparision', fontsize=20)

enhancersFig.savefig("images/enhancers.png") 
fig27.figure.clf()


# In[70]:


fig37 = finalConcat.groupby([9]).nunique()[7].plot(kind='bar',color="red",title='Genes', ylabel='Unique Genes Count',
         xlabel='Method', figsize=(10, 10))
fig37.figure.suptitle('Gene comparision', fontsize=20)

fig37.figure.savefig("images/GeneComparsion.png") 
fig37.figure.clf()


# In[71]:


ax1=finalConcat.groupby([9]).nunique()[3].plot(kind='bar' , color="indigo",title='Enhancers',
          figsize=(6, 5))
ax2=finalConcat.groupby([9]).nunique()[7].plot(kind='bar',color="red",title='Genes',
         figsize=(6, 5))
indigo =mpatches.Patch(color='indigo', label='Enhancer')

red =mpatches.Patch(color='red', label='Gene')

plt.legend(handles=[indigo, red])
plt.title("Enhancer gene bar")
# plt.show()
plt.clf()


# In[72]:


import matplotlib.pyplot as plt1
ax1=finalConcat.groupby([9]).nunique()[3].plot(kind='bar' , color="indigo",title='Enhancers',logy=True,
          figsize=(6, 5))
ax2=finalConcat.groupby([9]).nunique()[7].plot(kind='bar',color="red",title='Genes',logy=True,
         figsize=(6, 5))
indigo =mpatches.Patch(color='indigo', label='Enhancer')

red =mpatches.Patch(color='red', label='Gene')
plt1.legend(handles=[indigo, red])
plt1.savefig('images/enhancerGene.png')
plt1.clf()


# In[73]:


allDF = chiaPet.merge(eqTL,on=3).merge(distance,on=3)


# In[74]:


len(allDF)


# In[75]:


allDF.head()


# In[76]:


allDF1 = allDF[(allDF['7_x'] == allDF['7_y']) & (allDF['7_x'] ==allDF[7])]


# In[77]:


allDF1.head()


# In[78]:


k = allDF1[['8_x','8_y']]


# In[79]:


k = k.rename(columns={'8_x': 'chiaPet', '8_y': 'eqtl'})


# In[80]:


k.head()


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

# In[81]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn3


# In[82]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn3

set1 = set(chiaPet[7])
set2 = set(distance[7])
set3 = set(eqTL[7])
fig =plt.figure(figsize=(10,10))

venn3([set1, set2, set3], ('chiapet', 'distance', 'eqtl'))
fig.suptitle('Venn diagram of unique Genes', fontsize=20)
plt.savefig('images/AllGeneComparsion.png')
plt.clf()


# In[83]:


set1 = set(chiaPet[3])
set2 = set(distance[3])
set3 = set(eqTL[3])
fig = plt.figure(figsize=(10,10))
fig.suptitle('Venn diagram of number of unique Enhancer', fontsize=20)

venn3([set1, set2, set3], ('chiapet', 'distance', 'eqtl'))
plt.savefig('images/AllEnhancerComparsion.png')
plt.clf()


# Enhancer-Gene combination used

# In[84]:


from collections import Counter


# In[85]:


set1 = set(chiaPet[10])
set3 = set(distance[10])
set2 = set(eqTL[10])
sets = [set1,set2,set3]
fig = plt.figure(figsize=(10,10))
fig.suptitle('Venn diagram of number of unique Enhancer-Gene', fontsize=20)
ax = plt.gca()
v = venn3([set1, set2, set3], ('chiaPet', 'eqTL', 'distance'), ax = ax)
h, l = [],[]
for i in ['10','01','111']:
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

plt.savefig('images/enhancerGeneVenn.png')
plt.clf()


# In[86]:


# sns.set(rc={'figure.figsize':(11.7,8.27)})

# sns.lineplot(data=k)


# In[87]:


from matplotlib_venn import venn2


# In[88]:


set1 = set(chiaPet[10])
set2 = set(distance[10])
set3 = set(eqTL[10])
fig = plt.figure(figsize=(10,10))
fig.suptitle('Venn diagram of number of unique Enhancer', fontsize=20)

ax = venn2([set1,  set3], ('chiapet',  'eqtl'))
# plt.show()
plt.savefig('images/enhancerchiaPetEQTL.png')
plt.clf()


# In[89]:


chiaPet.head()


# In[90]:


from collections import Counter


# In[91]:


df1 = dict(Counter(chiaPet[3]))


# In[92]:


df2 = dict(Counter(distance[3]))


# In[93]:


df3 = dict(Counter(eqTL[3]))


# In[94]:


newdf =  pd.DataFrame.from_dict(df1, orient="index").reset_index()


# In[95]:


newdf=newdf.rename(columns={0: 'repetitions'})


# In[96]:


newdf[1]= 'chiaPet'


# In[97]:


newdf.head()


# In[ ]:





# In[98]:


newdf_distance =  pd.DataFrame.from_dict(df2, orient="index").reset_index()


# In[ ]:





# In[99]:


newdf_distance = newdf_distance.rename(columns={0: 'repetitions'})


# In[100]:


newdf_distance[1]= 'distance'


# In[101]:


newdf_distance.head()


# In[102]:


newdf_eqtl =  pd.DataFrame.from_dict(df3, orient="index").reset_index()


# In[103]:


newdf_eqtl = newdf_eqtl.rename(columns={0: 'repetitions'})


# In[104]:


newdf_eqtl[1]= 'eqtl'


# In[ ]:





# In[105]:


del(df1)
del(df2)
del(df3)
del(chiaPet)
del(distance)
del(eqTL)


# In[107]:


fig22 = sns.histplot(data=newdf, x='repetitions')

chiapetFig = fig22.get_figure()
chiapetFig.suptitle('Histogram for chiapet of repetitions of Enhancer-Gene', fontsize=20)

chiapetFig.savefig("images/chiapetHisto.png") 


# In[108]:


newdf_distance[newdf_distance.repetitions  > 1].count()


# In[109]:


newdf_distance.head()


# In[110]:


fig223=sns.histplot(data=newdf_distance, x='repetitions',bins=5)

distanceFig = fig223.get_figure()
distanceFig.suptitle('Histogram for Distance of repetitions of Enhancer-Gene', fontsize=20)

distanceFig.savefig("images/distanceHisto.png") 
fig223.figure.clf()


# In[112]:


fig228 = sns.histplot(data=newdf_eqtl, x='repetitions')
eqtlFig = fig228.get_figure()
eqtlFig.suptitle('Histogram for Distance of repetitions of Enhancer-Gene', fontsize=20)

eqtlFig.savefig("images/eqtlHisto.png") 
fig228.figure.clf()


# In[ ]:





# # P value analysis(common sets)
# # distance as qualitative in distance method
# # P value display on diagram
# # filter more data (remove dots)
