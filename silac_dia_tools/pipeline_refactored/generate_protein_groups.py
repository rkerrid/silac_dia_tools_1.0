# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:49:28 2024

@author: robbi
"""
import pandas as pd
import numpy as np


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
        # replace precursor quantity with summed silac channels as input for direct lefq and as 'total intensity' for href quantification
        merged_df['Precursor.Quantity'] = merged_df['Precursor.Translated H'] + merged_df['Precursor.Translated M'] + merged_df['Precursor.Translated L'] 
    else:
        print("Column 'Label H' not found in DataFrame")
        print(merged_df.columns.values.tolist())
    
    return merged_df

def split_data_by_intensity_type(df):
    # ms1 h/l and h/m as well as precursor h/l h/m (perhaps separate dfs??)
    precursor_translated_df = df.copy(deep=True)
    precursor_translated_df = precursor_translated_df[['Run','Protein.Group', 'Precursor.Id', 'Precursor.Quantity', 'Precursor.Translated H', 'Precursor.Translated M', 'Precursor.Translated L']]
    
    ms1_translated_df = df.copy(deep=True)
    ms1_translated_df = ms1_translated_df[['Run','Protein.Group', 'Precursor.Id', 'Precursor.Quantity', 'Ms1.Translated H', 'Ms1.Translated M', 'Ms1.Translated L']]
    
    return precursor_translated_df, ms1_translated_df

def calculate_precursor_ratios(df, quantification):
    
    df[f'{quantification} H/T'] = df[f'{quantification} H'] / df['Precursor.Quantity']
    df[f'{quantification} M/T'] = df[f'{quantification} M'] / df['Precursor.Quantity']
    df[f'{quantification} L/T'] = df[f'{quantification} L'] / df['Precursor.Quantity']
    df['Lib.PG.Q.Value'] = 0
    return df

def compute_protein_level(df): # this function looks for at least 3 valid values for each ratio and sums Precursor.Quantity (which is the sum of precursor translated values) for total intensity
    
    # Function to filter values and compute median
    def valid_median(series):
        valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
        valid_series = np.log2(valid_series)
        return 2 **valid_series.median() if len(valid_series) >= 2 else np.nan
    
    def valid_sum(series):
        valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
        return valid_series.sum() 
    
    result = df.groupby(['Protein.Group', 'Run']).agg({
        'Precursor.Translated H/T': valid_median,
        'Precursor.Translated M/T': valid_median,
        'Precursor.Translated L/T': valid_median,
        'Precursor.Quantity': valid_sum        
    })
    result['H'] = result['Precursor.Translated H/T']*result['Precursor.Quantity']
    result['M'] = result['Precursor.Translated M/T']*result['Precursor.Quantity']
    result['L'] = result['Precursor.Translated L/T']*result['Precursor.Quantity']
    result = result.reset_index()
    
    cols = ['Run', 'Protein.Group', 'Precursor.Quantity', 'H', 'M', 'L'] # fix this 
    return result[cols]


# for dlfq
def format_silac_channels_dlfq(df):
    # Pivot for each label
    pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' L')
    pivot_M = df[df['Label'] == 'M'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' M')
    
    # Merge the pivoted DataFrames
    merged_df = pd.concat([pivot_L, pivot_M], axis=1)
    # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
    merged_df.reset_index(inplace=True)

    # remove all rows that contain a NaN under the Label H column (i.e., no H precursor is present for that row)
    if 'Label H' in merged_df.columns:
        # Apply dropna on merged_df instead of df
        merged_df = merged_df.dropna(subset=['Precursor.Translated H'])
        # replace precursor quantity with summed silac channels as input for direct lefq and as 'total intensity' for href quantification
        merged_df['Precursor.Quantity'] = merged_df['Precursor.Translated H'] + merged_df['Precursor.Translated M'] + merged_df['Precursor.Translated L'] 
    else:
        print("Column 'Label H' not found in DataFrame")
        print(merged_df.columns.values.tolist())
        merged_df['Precursor.Quantity'] = merged_df['Precursor.Translated M'] + merged_df['Precursor.Translated L'] 
    return merged_df

def split_data_by_intensity_type_dlfq(df):
    # ms1 h/l and h/m as well as precursor h/l h/m (perhaps separate dfs??)
    precursor_translated_df = df.copy(deep=True)
    precursor_translated_df = precursor_translated_df[['Run','Protein.Group', 'Precursor.Id', 'Precursor.Quantity', 'Precursor.Translated M', 'Precursor.Translated L']]
    
    ms1_translated_df = df.copy(deep=True)
    ms1_translated_df = ms1_translated_df[['Run','Protein.Group', 'Precursor.Id', 'Precursor.Quantity', 'Ms1.Translated M', 'Ms1.Translated L']]
    
    return precursor_translated_df, ms1_translated_df

def calculate_precursor_ratios_dlfq(df, quantification):
    df[f'{quantification} M/T'] = df[f'{quantification} M'] / df['Precursor.Quantity']
    df[f'{quantification} L/T'] = df[f'{quantification} L'] / df['Precursor.Quantity']
    df['Lib.PG.Q.Value'] = 0
    return df

def compute_protein_level_dlfq(df): # this function looks for at least 3 valid values for each ratio and sums Precursor.Quantity (which is the sum of precursor translated values) for total intensity
    
    # Function to filter values and compute median
    def valid_median(series):
        valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
        valid_series = np.log2(valid_series)
        return 2 **valid_series.median() if len(valid_series) >= 2 else np.nan
    
    def valid_sum(series):
        valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
        return valid_series.sum() 
    
    result = df.groupby(['Protein.Group', 'Run']).agg({
        'Precursor.Translated M/T': valid_median,
        'Precursor.Translated L/T': valid_median,
        'Precursor.Quantity': valid_sum        
    })
    result['M'] = result['Precursor.Translated M/T']*result['Precursor.Quantity']
    result['L'] = result['Precursor.Translated L/T']*result['Precursor.Quantity']
    result = result.reset_index()
    
    cols = ['Run', 'Protein.Group', 'Precursor.Quantity', 'M', 'L', 'Precursor.Translated M/T', 'Precursor.Translated L/T' ] # fix this 
    return result[cols]



  