# -*- coding: utf-8 -*-
"""
Created on Sun Apr 18 11:32:44 2021

@author: Dell
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
import scipy.cluster.hierarchy as shc

from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
from sklearn.metrics import calinski_harabasz_score
from sklearn.metrics import davies_bouldin_score
from sklearn.decomposition import PCA

from munkres import Munkres

from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import plot_tree
from sklearn.tree import export_graphviz
from sklearn.ensemble import RandomForestClassifier
import graphviz

#returns dictionary of class labels
def get_class_labels(class_array):
    count = 0
    classes = {}
    for name in class_array:
        if name not in classes.keys():
            classes[name] = count
            count += 1   
    return classes
def set_labels_numerical(class_array, labels_dict):
    num_class_array = class_array.copy()
    for i in range(0,class_array.size):
        num_class_array[i] = labels_dict[class_array[i]]
    return num_class_array
def generate_confusion_matrix(classes, labels):
    labels = labels[:320]
    confusion_matrix = np.zeros((8,8))
    for name,label in zip(classes,labels):
        confusion_matrix[name][label] += 1   
    return confusion_matrix
def get_cost_matrix (matrix):
    cost_matrix = []
    for row in matrix:
        cost_row = []
        for col in row:
            cost_row += [1000000 - col]
        cost_matrix += [cost_row]
    return cost_matrix
def reverse_dict(dictionary):
    new_dict = {}
    for key in dictionary:
        new_dict[dictionary[key]] = key
    return new_dict
def get_cluster_label_dict(indexes,class_labels_dict):
    rev_class_label_dict = reverse_dict(class_labels_dict)
    cluster_label_dict = {}
    for row,column in indexes:
        cluster_label_dict[rev_class_label_dict[row]] = column
    return cluster_label_dict
def score_clustering(class_labels,kmeans_labels,method,unk_present = False):
    class_labels_dict = get_class_labels(class_labels)
    num_labels = set_labels_numerical(class_labels, class_labels_dict)
    confusion_matrix = generate_confusion_matrix(num_labels,kmeans_labels)
    cost_matrix = get_cost_matrix(confusion_matrix)
    m = Munkres()
    indexes = m.compute(cost_matrix)
    total = 0
    for row, column in indexes:
        value = confusion_matrix[row][column]
        total += value
    print(method,"clustering accurately identified:", total/len(kmeans_labels), "of protein locals")
    cluster_label_dict = get_cluster_label_dict(indexes,class_labels_dict)
    return reverse_dict(cluster_label_dict)
    
ftest = "project_ecoli_localization.test.csv"
ftrain = "project_ecoli_localization.train.csv"

df_train = pd.read_csv(ftrain)
df_test = pd.read_csv(ftest)
#must be this order for concat
full_df = df_train.append(df_test, ignore_index = True)
class_train = np.array(df_train['class'])
df_train = df_train.drop('class', axis = 1)
class_test = np.array(df_test['class'])
df_test = df_test.drop('class', axis = 1)
class_full = np.array(full_df['class'])
print(full_df)
full_df = full_df.drop('class', axis = 1)
# represent using pair plot
""" g = sns.pairplot(full_df, height=2.0, kind='reg')
plt.savefig('pairplot_proj2.png')
plt.show() """


#find optimal number of kmeans clusters
kmeans_wcss = []
silhouette_scores = []#maximize
CH_scores = []#maximize
DB_scores = []#minimize
for i in range(2,16):
    kmeans_model = KMeans(n_clusters = i,random_state=42)
    kmeans_model = kmeans_model.fit(df_train)
    kmeans_labels = np.array(kmeans_model.labels_)
    kmeans_wcss.append(kmeans_model.inertia_)
    silhouette_scores.append(silhouette_score(df_train,kmeans_labels))
    CH_scores.append(calinski_harabasz_score(df_train,kmeans_labels))
    DB_scores.append(davies_bouldin_score(df_train,kmeans_labels))

""" fig1,ax = plt.subplots(4,1,sharex = True, figsize=(10,10))
ax[0].scatter(range(2,16),silhouette_scores)
ax[0].set_ylabel("sil_score")
ax[1].scatter(range(2,16),CH_scores)
ax[1].set_ylabel("ch score")
ax[2].scatter(range(2,16),DB_scores)
ax[2].set_ylabel("db score")
ax[0].set_title('kmeans scores')
ax[3].scatter(range(2,16),kmeans_wcss)
ax[3].set_ylabel("wcss")
plt.show() """

kmeans_model = KMeans(n_clusters = 8,random_state=42)
kmeans_model = kmeans_model.fit(df_train)
kmeans_labels = np.array(kmeans_model.labels_)

score_clustering(class_train,kmeans_labels,"kmeans")

# fig4,ax4 = plt.subplots()
# ax4.scatter(df_train['mcg'], df_train['alm2'], c=kmeans_labels, cmap = 'hsv')
# ax4.set_xlabel('mcg')
# ax4.set_ylabel('alm2')
# ax4.set_title("cell loc kmeans", fontsize=14)
# plt.savefig('cell loc kmeans.png')
# plt.show()

#GMM
gmm = GaussianMixture(n_components=4, random_state=42)
gmm_labels = gmm.fit_predict(df_train)
gmm_labels = np.array(gmm_labels)
score_clustering(class_train,gmm_labels,"gmm")

#ahc
#find optimal number of ahc clusters
ahc_silhouette_scores = []#maximize
ahc_CH_scores = []#maximize
ahc_DB_scores = []#minimize
for i in range(2,16):
    ahc = AgglomerativeClustering(n_clusters = i, affinity='euclidean', linkage='complete')
    ahc = ahc.fit(df_train)
    ahc_labels = np.array(ahc.labels_)
    ahc_silhouette_scores.append(silhouette_score(df_train,ahc_labels))
    ahc_CH_scores.append(calinski_harabasz_score(df_train,ahc_labels))
    ahc_DB_scores.append(davies_bouldin_score(df_train,ahc_labels))

""" fig2,ax2 = plt.subplots(3,1,sharex = True, figsize=(10,10))
ax2[0].scatter(range(2,16),ahc_silhouette_scores)
ax2[0].set_ylabel("ahc_sil_score")
ax2[1].scatter(range(2,16),ahc_CH_scores)
ax2[1].set_ylabel("ahc_ch score")
ax2[2].scatter(range(2,16),ahc_DB_scores)
ax2[2].set_ylabel("ahc_db score")
ax2[0].set_title('ahc scores complete')
plt.show() """

#complete linkage labels
ahc_comp = AgglomerativeClustering(n_clusters = 8, affinity='euclidean', linkage='complete')
ahc_comp = ahc_comp.fit(df_train)
ahc_complete_labels = np.array(ahc_comp.labels_)
score_clustering(class_train,ahc_complete_labels,"AHC_complete")
#average linkage labels
ahc_avg = AgglomerativeClustering(n_clusters = 8, affinity='euclidean', linkage='average')
ahc_avg = ahc_avg.fit(df_train)
ahc_avg_labels = np.array(ahc_avg.labels_)
score_clustering(class_train,ahc_avg_labels,"AHC_avg")

# fig3,ax3 = plt.subplots()
# ax3.scatter(df_train['mcg'], df_train['alm2'], c=ahc_avg_labels, cmap = 'hsv')
# ax3.set_xlabel('mcg')
# ax3.set_ylabel('alm2')
# ax3.set_title("cell loc ahc avg", fontsize=14)
# plt.savefig('cell loc ahc avg.png')
# plt.show()
""" plt.figure(figsize=(10, 7))
plt.title("Dendogram for Protein Localization")
dend = shc.dendrogram(shc.linkage(df_train, method='complete'))
plt.savefig('dendogram_complete.png')
plt.show() """

#predict using AHC w/ avg linkage
ahc_avg = AgglomerativeClustering(n_clusters = 8, affinity='euclidean', linkage='average')
ahc_avg = ahc_avg.fit(full_df)
ahc_avg_labels = np.array(ahc_avg.labels_)
cluster_identities = score_clustering(class_train,ahc_avg_labels,"AHC_avg")
unk_labels = []
for unk in ahc_avg_labels[321:]:
    unk_labels.append(cluster_identities[unk])
print("unknown labels:",unk_labels)


#decision tree
""" tree = DecisionTreeClassifier()
tree = tree.fit(df_train,class_train)
tree_dot = export_graphviz(tree, out_file=None, 
                feature_names=full_df.columns,  
                class_names=list(cluster_identities),  
                filled=True, rounded=True,  
                special_characters=True)
graph = graphviz.Source(tree_dot)
graph.render("cell location_tree") """

#predict using random forest
forest = RandomForestClassifier(oob_score=True).fit(df_train,class_train)
print("random forest score = ", forest.oob_score_)
print("feature weights",forest.feature_importances_)
unk_labels = forest.predict(df_test)
print("random forest unknown labels:",unk_labels)

# PCA
pca = PCA(random_state = 42)
pca.fit(df_train)
# print(pca.components_)
pca_train = pca.transform(df_train)
kmeans_model = KMeans(n_clusters = 8,random_state=42)
kmeans_model = kmeans_model.fit(pca_train)
kmeans_labels = np.array(kmeans_model.labels_)

score_clustering(class_train,kmeans_labels,"kmeans with pca")
#complete linkage labels
ahc_comp = AgglomerativeClustering(n_clusters = 8, affinity='euclidean', linkage='complete')
ahc_comp = ahc_comp.fit(pca_train)
ahc_complete_labels = np.array(ahc_comp.labels_)
score_clustering(class_train,ahc_complete_labels,"AHC_complete with pca")
#average linkage labels
ahc_avg = AgglomerativeClustering(n_clusters = 8, affinity='euclidean', linkage='average')
ahc_avg = ahc_avg.fit(pca_train)
ahc_avg_labels = np.array(ahc_avg.labels_)
score_clustering(class_train,ahc_avg_labels,"AHC_avg with pca")