# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 18:00:04 2021

@author: Dell
"""

import os
import pandas as pd
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
