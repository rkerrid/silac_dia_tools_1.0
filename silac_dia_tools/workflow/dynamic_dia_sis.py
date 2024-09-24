# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 10:43:55 2024

@author: robbi
"""
import pandas as pd
import numpy as np
import time
from tqdm import tqdm
from icecream import ic
from .utils import manage_directories


class DynamicDiaSis:
    def __init__(self, path, filtered_report):
        self.path = path
        self.filtered_report = filtered_report
        self.update = True
        
        self.formatted_precursors = None
        self.protein_groups = None
        
    def generate_protein_groups(self):
        print('Begin protein groups')
        precursors = self.get_filtered_set(self.filtered_report)
        
        href = self.generate_href_df( precursors)
        LH_precursors, MH_precursors = self.get_precursor_ratios( precursors)
        
        L_pg = self.get_protein_ratios(LH_precursors, 'L')
        M_pg = self.get_protein_ratios(MH_precursors, 'M')
            
        protein_groups = self.merge_data(href, L_pg, M_pg)
        
        return protein_groups
        
    def get_filtered_set(self, df):
        df['filter_passed_L'] = df['filter_passed_L'].astype(bool)
        df['filter_passed_H'] = df['filter_passed_H'].astype(bool)
        df['filter_passed_M'] = df['filter_passed_M'].astype(bool)

        return  df[(df['filter_passed_H'])]
        
    def generate_href_df(self, precursors):
        df_href_precursors = precursors[['Run', 'protein_group', 'precursor_id', 'precursor_quantity_H']]
        df_href_precursors = df_href_precursors.replace(0, np.nan)
        df_href_precursors = df_href_precursors.dropna(subset=['precursor_quantity_H'])
        
        # Group by protein group and apply the custom aggregation
        grouped = df_href_precursors.groupby(['Run', 'protein_group']).agg({
                'precursor_quantity_H':'sum'
        }
        ).reset_index()
        
        href_df = grouped.groupby('protein_group').agg({
            'precursor_quantity_H':'median'
        })
        
        grouped = grouped.rename(columns={'precursor_quantity_H':'H'})
        href_df = href_df.rename(columns={'precursor_quantity_H':'H_ref'})
        
        href_df = pd.merge(grouped, href_df, on=['protein_group'])
        return href_df
    
    def get_precursor_ratios(self, precursors):
        precursors = precursors.copy()
        
        precursors['precursor_translated_LH'] = precursors['precursor_translated_L'] / precursors['precursor_translated_H'] 
        precursors['ms1_translated_LH'] = precursors['ms1_translated_L'] / precursors['ms1_translated_H'] 
        precursors['precursor_quantity_LH'] = precursors['precursor_quantity_L'] / precursors['precursor_quantity_H'] 

        
        precursors['precursor_translated_MH'] = precursors['precursor_translated_M'] / precursors['precursor_translated_H'] 
        precursors['ms1_translated_MH'] = precursors['ms1_translated_M'] / precursors['ms1_translated_H']
        precursors['precursor_quantity_MH'] = precursors['precursor_quantity_M'] / precursors['precursor_quantity_H'] 

        
        columns_to_replace = ['precursor_translated_LH', 'ms1_translated_LH', 'precursor_translated_MH', 'ms1_translated_MH', 'precursor_quantity_LH', 'precursor_quantity_MH']
        precursors[columns_to_replace] = precursors[columns_to_replace].replace([np.inf, -np.inf, 0.0], np.nan)
        
        index_cols = ['Run','protein_group','precursor_id']
        
        LH_precursors = precursors[index_cols + ['precursor_translated_LH', 'ms1_translated_LH', 'precursor_quantity_LH']]
        LH_precursors = LH_precursors.dropna(subset=['precursor_translated_LH', 'precursor_quantity_LH'])
        
        MH_precursors = precursors[index_cols + ['precursor_translated_MH', 'ms1_translated_MH','precursor_quantity_MH']]
        MH_precursors = MH_precursors.dropna(subset=['precursor_translated_MH', 'precursor_quantity_MH'])
        
        return LH_precursors, MH_precursors 

    def get_protein_ratios(self, df, channel):
        runs = df['Run'].unique()
        runs_list = []
                
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
            
        
            def combined_median(ms1_series, precursor_series, quantity_series):
                if len(quantity_series.dropna()) <= 1:  # Remove NaNs before counting
                    return np.nan
                else:
                    combined_series = np.concatenate([ms1_series, precursor_series, quantity_series])
                    combined_series = combined_series[~np.isnan(combined_series)]
                    combined_series = np.log2(combined_series)  # Log-transform the combined series
                    return 2**np.median(combined_series)  # Return the median of the log-transformed values
             
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby(['protein_group']).apply(lambda x: pd.Series({
                f'{channel}H_ratio': combined_median(x[f'ms1_translated_{channel}H'], 
                                                     x[f'precursor_translated_{channel}H'], 
                                                     x[f'precursor_quantity_{channel}H'])})).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
        
        result = pd.concat(runs_list, ignore_index=True)
        
        return result
    
    def merge_data(self, href, L_pg, M_pg):
        
        protein_ratios = pd.merge(L_pg, M_pg, on=['Run', 'protein_group'], how='outer' )
        protein_ratios = pd.merge(protein_ratios, href, on=['protein_group', 'Run'], how='left')
        protein_ratios['L'] = protein_ratios['LH_ratio'] * protein_ratios['H_ref']
        protein_ratios['M'] = protein_ratios['MH_ratio'] * protein_ratios['H_ref']
        
        return protein_ratios
   