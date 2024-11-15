# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 15:07:41 2024

@author: robbi
"""

import pandas as pd
import numpy as np
import time
from tqdm import tqdm
from icecream import ic
from .utils import manage_directories

from silac_dia_tools.workflow.utils import dlfq_functions as dlfq


class DynamicSilac:    
    def __init__(self, path, filtered_report):
        self.path = path
        self.filtered_report = filtered_report
        
        self.protein_groups = None
    
    def generate_protein_groups(self):  
        start_time = time.time()        
        self.filtered_report = self.filter_data(self.filtered_report)
        
        precursor_ratios = self.calculate_precursor_ratios(self.filtered_report)
        protein_group_ratios = self.compute_protein_level_ratios(precursor_ratios)
        protein_group_ratios.to_csv(f'{self.path}/preprocessing/protein_group_ratios.csv', sep=',')
        protein_intensities_dlfq = self.perform_lfq(precursor_ratios)
        # protein_intensites_unnormalized = self.get_unnormalized_intensities(precursor_ratios)
        
        self.protein_groups = self.merge_data(protein_group_ratios, protein_intensities_dlfq)
        
        self.protein_groups = self.extract_M_and_L(self.protein_groups)
        
        end_time = time.time()
        print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        return self.protein_groups
    
    def filter_data(self, df):
        df['filter_passed_L'] = df['filter_passed_L'].astype(bool)
        df['filter_passed_pulse'] = df['filter_passed_pulse'].astype(bool)
        return  df[(df['filter_passed_L']) | (df['filter_passed_pulse'])]
    
    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        df = df.copy()
        # cols = ['precursor_quantity_L', 'precursor_translated_L', 'ms1_translated_L', 'precursor_translated_pulse', 'ms1_translated_pulse', 'precursor_quantity_L', 'precursor_quantity_pulse']
        
        df.loc[:, 'Precursor.Quantity'] = df['precursor_quantity_L'].fillna(0) + df['precursor_quantity_pulse'].fillna(0)
        
        df.loc[:, 'precursor_quantity_pulse_L_ratio'] = df['precursor_quantity_pulse'] / df['precursor_quantity_L'] 
        df.loc[:, 'precursor_translated_pulse_L_ratio'] = df['precursor_translated_pulse'] / df['precursor_translated_L'] 
        df.loc[:, 'ms1_translated_pulse_L_ratio'] = df['ms1_translated_pulse'] / df['ms1_translated_L']
        
        df.loc[:, 'precursor_quantity_pulse_L_ratio'] = df['precursor_quantity_pulse_L_ratio'].replace([np.inf, 0], np.nan)
        df.loc[:, 'precursor_translated_pulse_L_ratio'] = df['precursor_translated_pulse_L_ratio'].replace([np.inf, 0], np.nan)
        df.loc[:, 'ms1_translated_pulse_L_ratio'] = df['ms1_translated_pulse_L_ratio'].replace([np.inf, 0], np.nan)
        
        df.loc[:, 'Lib.PG.Q.Value'] = 0
        
        return df
    
    def compute_protein_level_ratios(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
    
        df = df.dropna(subset=['precursor_translated_pulse_L_ratio','precursor_quantity_pulse_L_ratio'])
    
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
            
            def combined_median(pre_translated, pre_quantity, ms1_translated):
                
                if len(pre_quantity.dropna()) <= 1:  # Remove NaNs before counting
                    return np.nan
                else:
                    combined_series = np.concatenate([pre_translated, pre_quantity, ms1_translated])
                    combined_series = combined_series[~np.isnan(combined_series)]
                    combined_series = np.log2(combined_series)  # Log-transform the combined series
                    return 2**np.median(combined_series)  # Return the median of the log-transformed values
    
            # Group by protein group and apply the custom aggregation
            grouped_run = run_df.groupby(['protein_group']).apply(lambda x: pd.Series({
                'pulse_L_ratio': combined_median(x['precursor_translated_pulse_L_ratio'], x['precursor_quantity_pulse_L_ratio'], x['ms1_translated_pulse_L_ratio'])
            })).reset_index()
    
            grouped_run['Run'] = run
            runs_list.append(grouped_run)
    
        result = pd.concat(runs_list, ignore_index=True)
        return result

    def perform_lfq(self, df):
            manage_directories.create_directory(self.path, 'directLFQ_output')
            df = df[['Run', 'protein_group', 'precursor_id', 'Precursor.Quantity', 'Lib.PG.Q.Value']]
            df = df.rename(columns={'protein_group':'Protein.Group', 'precursor_id':'Precursor.Id'})
            path = f'{self.path}directLFQ_output/'
            df.to_csv(f'{path}dflq_formatted_report.tsv', sep='\t')
            dlfq_output_file = f'{path}dlfq_protein_intensities.tsv'
            
            dlfq.run_lfq(f'{path}dflq_formatted_report.tsv', file=dlfq_output_file, num_cores=1)
            dlfq_df = pd.read_csv(dlfq_output_file, sep='\t')
           
            # Drop the 'Unnamed: 0' column
            dlfq_df = dlfq_df.drop(columns=['Unnamed: 0', 'protein'])
           
            # Melt the DataFrame
            result = pd.melt(dlfq_df, id_vars=['Protein.Group'], var_name='Run', value_name='Intensity')
            result = result.rename(columns={'Protein.Group':'protein_group','Intensity':'normalized_intensity'})
            result = result[result['normalized_intensity'] != 0]
            return result
    
    def get_unnormalized_intensities(self, df):
        df = df.groupby(['Run','protein_group']).agg({
            'Precursor.Quantity':'median'
        }).reset_index()
        
        return df
    
    def merge_data(self, protein_group_ratios, protein_intensities_dlfq):
        # protein_groups = pd.merge(protein_group_ratios, protein_intensites_unnormalized,on=['protein_group','Run'], how='left')
        protein_groups = pd.merge(protein_group_ratios, protein_intensities_dlfq,on=['protein_group','Run'], how='left')
        protein_groups = protein_groups.dropna(subset=['normalized_intensity'])
        return protein_groups
    
    def extract_M_and_L(self, df):
        
        def calculate_M_and_L(row):
           ratio = row['pulse_L_ratio']
           normalized_intensity = row['normalized_intensity']
        
           L_norm = normalized_intensity / (ratio + 1)
           pulse_norm = normalized_intensity - L_norm
            
           return pd.Series([pulse_norm, L_norm], index=['pulse', 'L'])
        
         
        df[['pulse', 'L']] = df.apply(calculate_M_and_L, axis=1)
        return df
    

    