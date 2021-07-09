# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 19:55:48 2021

@author: Dell
"""
import matplotlib.pyplot as plt
import seaborn as sns; sns.set() #set the style
import scipy.stats as stats
import numpy as np
import pandas as pd

from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from scipy import stats
from scipy.stats import ttest_ind

import statsmodels.api as sm
from statsmodels.graphics.gofplots import ProbPlot

filename = 'qm9.csv'
qm9 = pd.read_csv(filename, nrows = 150)
qm9 = qm9.drop(labels = ['mol_id','A', 'B', 'C', 'u0', 'h298',
'g298', 'u0_atom', 'h298_atom', 'g298_atom'], axis = 1)
print(qm9.loc[:,'mu':])

#pairplot
# g = sns.pairplot(qm9[1:], vars = qm9.columns[1:])
# g.map_upper(sns.regplot)
# g.map_lower(sns.kdeplot)
# g.savefig('pairplot.png')

#scale
scaled_labels = []
for i in range(len(qm9.columns[1:])):
    scaled_labels.append('scaled ' + qm9.columns[1:][i])
print(scaled_labels)    
qm9_scaled = pd.DataFrame(StandardScaler().fit_transform(qm9.loc[:,'mu':]),columns = scaled_labels)
qm9_scaled = pd.concat([qm9['smiles'],qm9_scaled],axis = 1)
#write to csv
qm9_to_csv = pd.concat([qm9,qm9_scaled],axis=1)
qm9_to_csv.to_csv('a.csv')

#PCA
#figs
fig1,ax1 = plt.subplots()
ax1.set_title('cum expl var ')
ax1.set_ylabel('expl var/var')
ax1.set_xlabel('# PCA')
fig3,ax3 = plt.subplots()
ax3.set_title('biplot ')
ax3.set_ylabel('PCA2')
ax3.set_xlabel('PCA1')
    
pca = PCA()
pca_fit = pca.fit(qm9_scaled.loc[:,'scaled mu':])
cumulative_exp_var_ratio = [0]
for i in range(len(pca_fit.explained_variance_ratio_)):
    cumulative_exp_var_ratio.append(pca_fit.explained_variance_ratio_[i] + cumulative_exp_var_ratio[i])

print(cumulative_exp_var_ratio)
ax1.plot(range(len(pca_fit.explained_variance_ratio_)+1),cumulative_exp_var_ratio)
plt.figure(fig1.number)
plt.savefig("cum expl var.png")

#biplot
PCs = ['PC1', 'PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10']
loadings = pd.DataFrame(pca.components_, columns=PCs, index=scaled_labels)
print(loadings)
#PCA

qm9_pca_scores = pd.DataFrame(pca_fit.transform(qm9_scaled.loc[:,'scaled mu':]),columns = PCs)
ax3.scatter(qm9_pca_scores['PC1'],qm9_pca_scores['PC2'])
plt.figure(fig3.number)
plt.savefig("biplot.png")

#part h
out_name = 'output.txt'
fout = open(out_name,'w+')
reg_targets = ['mu', 'cv', 'u298' , 'gap']
for c1,pc in enumerate(PCs[0:2]):
    for c2,t in enumerate(reg_targets):
        qm9_reg_fit = LinearRegression().fit([qm9_pca_scores[pc]],[qm9[t]])
        qm9_reg_predict = pd.DataFrame(qm9_reg_fit.predict([qm9_pca_scores[pc]]))
        print(r2_score(qm9[t],qm9_reg_predict.pivot()))
        # fout.write("r2_score: " + str(r2_score(qm9[t],qm9_reg_predict)))
        # fout.write("mse:" + str(mean_squared_error(r2_score(qm9[t],qm9_reg_predict))))
fout.close()
