#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
import os

file_path = sys.argv[1]
x = file_path[-22:-4] 

######## Calculation starts here ###############################
origin_data = pd.read_table(file_path, delimiter=",", index_col=0).fillna(0) #read a PCN table and fill nan with zero

                # Protein content network
PCN = origin_data # PCN without any normalization
PCN = PCN[PCN.columns[PCN.sum()>0]] # remove genus column if sum to zero
PCN_taxa_norm = PCN.apply(lambda x: x/sum(x), axis = 0).fillna(0) # Proteomic content network / normalized to each taxon
PCN_taxa_norm_T = PCN_taxa_norm.T
if_detected = PCN_taxa_norm_T > 0  # turn PCN numbers to 0 or 1
PCN_taxa_norm_T['num_detected'] = if_detected.sum(axis=1)  # obtain num_detected,
PCN_taxa_norm_all = PCN_taxa_norm_T.fillna(0).drop(columns='num_detected').T

                # create an empty dataframe
dij_matrix_norm = pd.DataFrame()
dij_list_norm = pd.DataFrame(
    index=range((len(PCN_taxa_norm_all.columns)) * (len(PCN_taxa_norm_all.columns) - 1) // 2), columns=['dij_list'])
dij_list_name_norm = pd.DataFrame(
    index=range((len(PCN_taxa_norm_all.columns)) * (len(PCN_taxa_norm_all.columns) - 1) // 2), columns=['VS'])
dij_num_norm = 0
                # getting a dij between any two bacteria genera
for i in range(0, len(PCN_taxa_norm_all.columns)):
    Gi = PCN_taxa_norm_all.iloc[:, i]
    print(i)
    for j in range(0, len(PCN_taxa_norm_all.columns)):
        Gj = PCN_taxa_norm_all.iloc[:, j]
        Gij = pd.concat([Gi, Gj], axis=1)
        dij = 1 - sum(Gij.apply(lambda x: min(x), axis=1)) / sum(Gij.apply(lambda x: max(x), axis=1))
        if i > j:
            dij_matrix_norm.at[Gij.columns.values[0], Gij.columns.values[1]] = dij
            dij_list_norm.at[dij_num_norm, 'dij_list'] = dij
            dij_num_norm = dij_num_norm + 1
            dij_list_name_norm.at[dij_num_norm-1] = Gij.columns.values[0] + '_vs_' + Gij.columns.values[1]
            dij_list_summary_norm = pd.concat([dij_list_name_norm, dij_list_norm], axis=1)
            dij_matrix_norm.to_csv('/beegfs/work/workspace/ws/ho_kezau83-snakemake-0/2_data_processing/3_PCN_Functional_distance/all/norm/matrix/all_norm_dij_matrix_'+ x+'.csv')
            dij_list_summary_norm.to_csv('/beegfs/work/workspace/ws/ho_kezau83-snakemake-0/2_data_processing/3_PCN_Functional_distance/all/norm/list/all_norm_dij_list_'+ x+'.csv')
