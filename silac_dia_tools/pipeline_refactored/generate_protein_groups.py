# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:49:28 2024

@author: robbi
"""
import pandas as pd

def format_silac_channels(df):
    # Pivot for each label
    pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' L')
    pivot_M = df[df['Label'] == 'M'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' M')
    pivot_H = df[df['Label'] == 'H'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' H')
    
    # Merge the pivoted DataFrames
    merged_df = pd.concat([pivot_L, pivot_M, pivot_H], axis=1)
    # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
    merged_df.reset_index(inplace=True)

    # remove all rows that contain a NaN under the Label H column (i.e., no H precursor is present for that row)
    if 'Label H' in merged_df.columns:
        # Apply dropna on merged_df instead of df
        merged_df = merged_df.dropna(subset=['Precursor.Translated H'])
    else:
        print("Column 'Label H' not found in DataFrame")
        print(merged_df.columns.values.tolist())
    
    return merged_df