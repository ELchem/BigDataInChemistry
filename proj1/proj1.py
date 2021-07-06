# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 22:28:59 2021

@author: Dell
"""

import sys


original_stdout = sys.stdout # Save a reference to the original standard output

import matplotlib.pyplot as plt
import seaborn as sns; sns.set() #set the style
import scipy.stats as stats
import numpy as np
import pandas as pd


def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

from sklearn import linear_model
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate

from scipy import stats
from scipy.stats import ttest_ind

import statsmodels.api as sm
from statsmodels.graphics.gofplots import ProbPlot

import os

#compound codes from paper
signaling = [681239, 
732517,
718781,
733504, 
715055, 
743414, 
745750, 
747599, 
747971, 
750690, 
683864]

DNA_Damaging = [3053,
409962,
138783,
125066,
750 ,
241240, 
79037,
3088 ,
119875,
26271 ,
45388 ,
13875 ,
109724 ,
8806 ,
45923 ,
26980 ,
24559 ,
762 ,
266046, 
25154 ,
77213 ,
14229 ,
85998 ,
362856 ,
6396 ,
9706 ,
34462]


Tubulin_directed = [628503, 
747973,
125973,
49842,
67574,
608210]

anthracyclines = [82151,
123127, 
256942 ,
141540,
256439 ,
616348 ,
301739 ,
122819,
609699 ,
246131]

antimetabolites = [102816, 
19893, 
755 ,
1390 ,
3590 ,
712807,
105014 ,
606869 ,
63878 ,
127716 ,
27640 ,
312887 ,
613327 ,
32065 ,
740 ,
686673, 
698037 ,
218321 ,
752] 

hormonal = [719344, 
23759 ,
88536,
12198 ,
702294 ,
10973 ,
713563 ,
719276 ,
719345 ,
71423 ,
38721,
23162,
747974, 
180973 ,
613680]

other = [296961,
18509,
706363,
719627,
169780,
369100,
747972,
177023,
113891,
747167,
630176,
66847,
122758,
701852,
721517]

def merge_list(df_list):
    df_merged = df_list[0]
    for df in df_list[1:]:
        df_merged = df_merged.merge(df,on = ['CellPanelName','CellLineName'])
    df_merged = df_merged.drop(['CellPanelName','CellLineName'],axis = 1)
    return df_merged

#stores DF by group
signaling_group = []
DNA_Damaging_group = []
hormonal_group = []
Tubulin_directed_group = []
other_group = []
anthracyclines_group = []
antimetabolites_group = []

#stores codes for each group
signaling_group_labels = []
DNA_Damaging_group_labels = []
hormonal_group_labels = []
Tubulin_directed_group_labels = []
other_group_labels = []
anthracyclines_group_labels = []
antimetabolites_group_labels = []

#list of label sets
label_list = [signaling_group_labels,DNA_Damaging_group_labels,hormonal_group_labels,
              Tubulin_directed_group_labels,other_group_labels,anthracyclines_group_labels,
              antimetabolites_group_labels]

with open('filename.txt', 'w') as f:
    sys.stdout = f # Change the standard output to the file we created.
    #reads files and assigns to groups
    #only reads label and logvalues
    #stores labels
    os.chdir('C:\\Users\\Dell\\Desktop\\Python Scripts\\big data in chem\\proj1\\data')
    for filename in os.listdir(os.getcwd()):
        data = pd.read_csv(filename, usecols = ['NSC','CellPanelName','CellLineName','logValue'])
        if int(data['NSC'][1]) == 715055:
            unknown = data
            data = data.rename(columns = {'logValue':'unknown'})
            signaling_group.append(data[['CellPanelName','CellLineName','unknown']])
            DNA_Damaging_group.append(data[['CellPanelName','CellLineName','unknown']])
            hormonal_group.append(data[['CellPanelName','CellLineName','unknown']])
            other_group.append(data[['CellPanelName','CellLineName','unknown']])
            anthracyclines_group.append(data[['CellPanelName','CellLineName','unknown']])
            antimetabolites_group.append(data[['CellPanelName','CellLineName','unknown']])
            Tubulin_directed_group.append(data[['CellPanelName','CellLineName','unknown']])
        elif int(data['NSC'][1]) in signaling:
            colname = str(data['NSC'][1])
            data =data.rename(columns = {'logValue':colname})
            signaling_group.append(data[['CellPanelName','CellLineName',colname]])
            signaling_group_labels.append(str(data['NSC'][1]))
            print(data.columns)
        elif int(data['NSC'][1]) in DNA_Damaging:
            colname = str(data['NSC'][1])
            data =data.rename(columns = {'logValue':colname})
            DNA_Damaging_group.append(data[['CellPanelName','CellLineName',colname]])
            DNA_Damaging_group_labels.append(str(data['NSC'][1]))
            print(data.columns)
        elif int(data['NSC'][1]) in hormonal:
            colname = str(data['NSC'][1])
            data =data.rename(columns = {'logValue':colname})
            hormonal_group.append(data[['CellPanelName','CellLineName',colname]])
            hormonal_group_labels.append(str(data['NSC'][1]))
            print(data.columns)
        elif int(data['NSC'][1]) in other:        
            colname = str(data['NSC'][1])
            data =data.rename(columns = {'logValue':colname})
            other_group.append(data[['CellPanelName','CellLineName',colname]])
            other_group_labels.append(str(data['NSC'][1]))
            print(data.columns)
        elif int(data['NSC'][1]) in anthracyclines:
            colname = str(data['NSC'][1])
            data =data.rename(columns = {'logValue':colname})
            anthracyclines_group.append(data[['CellPanelName','CellLineName',colname]])
            anthracyclines_group_labels.append(str(data['NSC'][1]))
            print(data.columns)
        elif int(data['NSC'][1]) in antimetabolites:
            colname = str(data['NSC'][1])
            data =data.rename(columns = {'logValue':colname})
            antimetabolites_group.append(data[['CellPanelName','CellLineName',colname]])
            antimetabolites_group_labels.append(str(data['NSC'][1]))
            print(data.columns)
        elif int(data['NSC'][1]) in Tubulin_directed:
            colname = str(data['NSC'][1])
            data =data.rename(columns = {'logValue':colname})
            Tubulin_directed_group.append(data[['CellPanelName','CellLineName',colname]])
            print(data.columns)
            Tubulin_directed_group_labels.append(str(data['NSC'][1]))
        else:
            print(data['NSC'][1],"not found")
    
    #forms group dataframes     
    sig_df = merge_list(signaling_group)
    dnaDam_df = merge_list(DNA_Damaging_group)
    horm_df = merge_list(hormonal_group)
    tub_df = merge_list(Tubulin_directed_group)
    other_df = merge_list(other_group)
    anth_df = merge_list(anthracyclines_group)
    anti_df = merge_list(antimetabolites_group)

    
    #list containing all datasets
    dataset_list = [sig_df,dnaDam_df,horm_df,tub_df,other_df,anth_df,anti_df]
    
    #remove rows with missing data
    for df in dataset_list:
        df.dropna(inplace = True)
    
        
    os.chdir('C:\\Users\\Dell\\Desktop\\Python Scripts\\big data in chem\\proj1')
    for count,df in enumerate(dataset_list):
        g = sns.pairplot(df)
        g.map_upper(sns.regplot)
        g.map_lower(sns.kdeplot)
        g.savefig("pairplot"+str(count)+".png")
        
    #list storing scipy linear regression models    
    LR_model_list = []
    OLS_model_list = []
    for count,df in enumerate(dataset_list):
        #OLS model
        modelOLS = sm.OLS(df['unknown'], sm.add_constant(df[label_list[count]]))
        model_fitOLS = modelOLS.fit()
        OLS_model_list.append(model_fitOLS)
        print(model_fitOLS.summary())
        for r in range(0,5):
            print("########\nrandom seed:",r,"\n########")
            #sklearn model
            training_set_x, test_set_x, training_set_y, test_set_y = train_test_split(df[label_list[count]],
                                                                                      df['unknown'], 
                                                                                      train_size=0.8,test_size = 0.2,
                                                                                      random_state=r)
            model_LR = LinearRegression(fit_intercept=True)
            model_LR_fit = model_LR.fit(training_set_x, training_set_y)
            # obtain R^2 for the training set
            print("Training set R^2= ",model_LR.score(training_set_x,training_set_y))
            
            # calculate RSE for the training set
            y_predict_train=model_LR.predict(training_set_x)
            RSS_train = np.sum((training_set_y - y_predict_train)**2)
            print("Training set MSE =",RSS_train/len(training_set_x))
            
            # obtain R^2 for the test set
            print("Test set R^2= ",model_LR.score(test_set_x,test_set_y))
            
            # calculate RSE for the test set
            y_predict_test=model_LR.predict(test_set_x)
            RSS_test = np.sum((test_set_y - y_predict_test)**2)
            print("Test set MSE =",RSS_train/len(training_set_x))
            
        #use cross validate method
        print("Using Cross Validate")
        score = cross_validate(model_LR,df[label_list[count]],df['unknown'],cv = 5, scoring = 'r2')
        print("metric = R^2",score)
        score = cross_validate(model_LR,df[label_list[count]],df['unknown'],cv = 5,scoring = 'neg_mean_squared_error')
        print("metric = MSE", score)
    sys.stdout = original_stdout 

















