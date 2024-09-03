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
        
        protein_intensities_dlfq = self.perform_lfq(precursor_ratios)
        protein_intensites_unnormalized = self.get_unnormalized_intensities(precursor_ratios)
        
        self.protein_groups = self.merge_data(protein_group_ratios, protein_intensites_unnormalized, protein_intensities_dlfq)
        
        self.protein_groups = self.extract_M_and_L(self.protein_groups)
        
        end_time = time.time()
        print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        return self.protein_groups
    
    def filter_data(self, df):
        return  df[(df['filter_passed_L']) | (df['filter_passed_pulse'])]
    
    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
     
        df.loc[:, 'Precursor.Quantity'] = df['precursor_quantity_L'].fillna(0) + df['precursor_quantity_pulse'].fillna(0)
        
        df.loc[:, 'precursor_translated_pulse/L'] = df['precursor_translated_pulse'] / df['precursor_translated_L'] 
        df.loc[:, 'ms1_translated_pulse/L'] = df['ms1_translated_pulse'] / df['ms1_translated_L']
        
        df.loc[:, 'precursor_translated_pulse/L'] = df['precursor_translated_pulse/L'].replace([np.inf, -np.inf], np.nan)
        df.loc[:, 'ms1_translated_pulse/L'] = df['ms1_translated_pulse/L'].replace([np.inf, -np.inf], np.nan)
        
        df.loc[:, 'Lib.PG.Q.Value'] = 0
        return df
    
    def compute_protein_level_ratios(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
    
        df = df.dropna(subset=['precursor_translated_pulse/L','ms1_translated_pulse/L'])
    
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
            
            def combined_median(ms1_series, precursor_series):
                # Combine the series
                combined_series = np.concatenate([ms1_series, precursor_series])
                
                # Check if the combined series contains only zeros
                if np.all(combined_series == 0):
                    return 0
                    
                quantifiable_ratios = combined_series[combined_series != 0] # at this point could select amount of ratios required for a valid ratio
                
                # Log-transform the modified series
                logged_ratios = np.log2(quantifiable_ratios)
                
                # Calculate the median in log space
                median_log = np.median(logged_ratios)
                
                # Back-transform from log space
                median_back_transformed = 2 ** median_log
                
                return median_back_transformed
    
            # Group by protein group and apply the custom aggregation
            grouped_run = run_df.groupby(['protein_group']).apply(lambda x: pd.Series({
                'pulse/L': combined_median(x['ms1_translated_pulse/L'], x['precursor_translated_pulse/L'])
            })).reset_index()
    
            grouped_run['Run'] = run
            runs_list.append(grouped_run)
    
        result = pd.concat(runs_list, ignore_index=True)
        return result

    def perform_lfq(self, df):
            df = df[['Run', 'protein_group', 'precursor_id', 'Precursor.Quantity', 'Lib.PG.Q.Value']]
            df = df.rename(columns={'protein_group':'Protein.Group', 'precursor_id':'Precursor.Id'})
            path = f'{self.path}'
            df.to_csv(f'{path}dflq_formatted_report.tsv', sep='\t')
            dlfq_output_file = f'{path}dlfq_protein_intensities.tsv'
            
            dlfq.run_lfq(f'{path}dflq_formatted_report.tsv', file=dlfq_output_file, num_cores=1)
            dlfq_df = pd.read_csv(dlfq_output_file, sep='\t')
           
            # Drop the 'Unnamed: 0' column
            dlfq_df = dlfq_df.drop(columns=['Unnamed: 0', 'protein'])
           
            # Melt the DataFrame
            result = pd.melt(dlfq_df, id_vars=['Protein.Group'], var_name='Run', value_name='Intensity')
            result = result.rename(columns={'Protein.Group':'protein_group','Intensity':'normalized_intensity'})
            return result
    
    def get_unnormalized_intensities(self, df):
        df = df.groupby(['Run','protein_group']).agg({
            'Precursor.Quantity':'median'
        }).reset_index()
        
        return df
    
    def merge_data(self, protein_group_ratios, protein_intensites_unnormalized, protein_intensities_dlfq):
        protein_groups = pd.merge(protein_group_ratios, protein_intensites_unnormalized,on=['protein_group','Run'], how='left')
        protein_groups = pd.merge(protein_groups, protein_intensities_dlfq,on=['protein_group','Run'], how='left')

        return protein_groups
    
    def extract_M_and_L(self, df):
        
        def calculate_M_and_L(row):
           ratio = row['pulse/L']
           normalized_intensity = row['normalized_intensity']
           intensity = row['Precursor.Quantity']
        
           L_norm = normalized_intensity / (ratio + 1)
           pulse_norm = normalized_intensity - L_norm
           
           L = intensity / (ratio + 1)
           pulse = intensity - L
            
            
           return pd.Series([pulse_norm, L_norm, pulse, L], index=['pulse_norm', 'L_norm', 'pulse', 'L'])
        
         
        df[['pulse_norm', 'L_norm', 'pulse', 'L']] = df.apply(calculate_M_and_L, axis=1)
        return df
    

    