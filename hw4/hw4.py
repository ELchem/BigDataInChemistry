# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 16:09:09 2021

@author: Dell
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score
from sklearn.metrics import davies_bouldin_score

def virus_to_number(virus):
    names = {'Hordevirus':1,'Tobravirus':2,'Tobamovirus':3,'Furovirus':4}
    return names[virus]
def target_to_df(target):
    nums_list = []
    names = {'Hordevirus':1,'Tobravirus':2,'Tobamovirus':3,'Furovirus':4}
    for i in target:
        nums_list.append(names[i])
    return pd.DataFrame(nums_list,columns = ['virus_number'])  
def number_to_virsu(numbers):
    names = {1:'Hordevirus', 2:'Tobravirus', 3:'Tobamovirus', 4:'Furovirus'}
    r = []
    for i in numbers:
        r.append(names[int(i)])
    return pd.DataFrame(r)
def assess_clustering(df):
    counts = df.value_counts(subset = ['virus_type'])
    #dictionary for storing virus, label pairs
    virus_labels = {}
    for v in counts.index:
        #v is a tuple
        virus = v[0] 
        #ordered count of labels 
        sub_counts = df[df['virus_type'] == virus].value_counts(subset = ['cluster'])
        #choose first label (higest count)
        cluster_label = sub_counts.keys()[0][0]
        
        virus_labels[virus] = cluster_label
    correct_count = 0
    for row in df.values:
        if virus_labels[row[0]] == row[1]:
            correct_count += 1
        else:
            continue
    print(virus_labels)
    return correct_count/df.shape[0]
        
    
FILE = 'prnn_viruses.csv'
df = pd.read_csv(FILE)

target = df['virus_type']
target_df = pd.DataFrame(df['virus_type'])
print(target)
df = df.drop(labels = 'virus_type',axis = 1)
print(df.mean(axis = 0))

# represent using pair plot
# g = sns.pairplot(df, height=2.0, kind='reg')
# plt.savefig('pairplot_shopping.png')
# plt.show()
fig1,ax1 = plt.subplots()
scatter = ax1.scatter(df['col_1'],df['col_2'],c=list(map(virus_to_number,target))
                      ,cmap = 'cool')
ax1.set_title('col_1 vs col_2')
ax1.set_ylabel('col_2')
ax1.set_xlabel('col_1')

#########KMEANS############
# data1 = df[['col_1','col_2']]
# df = df[['col_1','col_2']]
kmeans_score = []
silhouette_scores = []#maximize
CH_scores = []#maximize
DB_scores = []#minimize
for i in range(2,16):
    modelA = KMeans(n_clusters = i,random_state=42)
    model = modelA.fit(df)
    labels = np.array(model.labels_)
    kmeans_score.append(modelA.score(df))
    silhouette_scores.append(silhouette_score(df,labels))
    CH_scores.append(calinski_harabasz_score(df,labels))
    DB_scores.append(davies_bouldin_score(df,labels))

fig2,ax = plt.subplots(5,1,sharex = True, figsize=(10,10))
ax[0].scatter(range(2,16),silhouette_scores)
ax[0].set_ylabel("sil_score")
ax[1].scatter(range(2,16),CH_scores)
ax[1].set_ylabel("ch score")
ax[2].scatter(range(2,16),DB_scores)
ax[2].set_ylabel("db score")
ax[0].set_title('scores')
ax[3].scatter(range(2,16),kmeans_score)
ax[3].set_ylabel("kmeans score")

model = KMeans(n_clusters = 4, random_state=42).fit(df)
# print(model.labels_)
kmeans_model_labels = pd.DataFrame(model.labels_, columns=['cluster'])

virus_labels_kmeans = pd.concat([target_df,kmeans_model_labels], axis = 1)
print(virus_labels_kmeans)
print('score kmeans',assess_clustering(virus_labels_kmeans))
# virus_labels_kmeans.to_csv('kmeans_labels.csv')
##############GMM##############
gmm_model = GaussianMixture(4).fit_predict(df)
gmm_model_labels = pd.DataFrame(gmm_model, columns=['cluster'])
virus_labels_gmm = pd.concat([target_df,gmm_model_labels], axis = 1)
virus_labels_gmm.to_csv('gmm_cluster_labels.csv')
print('score gmm',assess_clustering(virus_labels_gmm))
#############agglomerative hier. clustering linkage = complete##########
agg_model = AgglomerativeClustering(4,linkage = 'complete').fit(df)
agg_model_labels = pd.DataFrame(agg_model.labels_, columns=['cluster'])
virus_labels_agg = pd.concat([target_df,gmm_model_labels], axis = 1)
virus_labels_agg.to_csv('agg_cluster_labels.csv')
print('score complete linkage',assess_clustering(virus_labels_agg))
#############agglomerative hier. clustering linkage = avg##########
agg_model = AgglomerativeClustering(4,linkage = 'average').fit(df)
agg_model_labels = pd.DataFrame(agg_model.labels_, columns=['cluster'])
virus_labels_agg = pd.concat([target_df,gmm_model_labels], axis = 1)
virus_labels_agg.to_csv('agg_cluster_labels_avg.csv')
print('score avg_linkage',assess_clustering(virus_labels_agg))