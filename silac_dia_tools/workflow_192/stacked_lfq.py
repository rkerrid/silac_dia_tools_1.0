# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:35:38 2024

@author: rkerrid
"""

import pandas as pd
import numpy as np
import time
from tqdm import tqdm
from icecream import ic
from .utils import manage_directories

from silac_dia_tools.workflow_192.utils import dlfq_functions as dlfq


class StackedLFQ:    
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
        protein_channel_mask = self.get_protein_level_channel_mask(precursor_ratios)
        
        # protein_intensites_unnormalized = self.get_unnormalized_intensities(precursor_ratios)
        ic(protein_group_ratios)
        ic(protein_channel_mask)
        ic(protein_intensities_dlfq)
        self.protein_groups = self.merge_data(protein_group_ratios, protein_channel_mask, protein_intensities_dlfq)
        ic(self.protein_groups)
        self.protein_groups = self.extract_M_and_L(self.protein_groups)
        
        end_time = time.time()
        print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        return self.protein_groups
    
    def filter_data(self, df):
        df['filter_passed_L'] = df['filter_passed_L'].astype(bool)
        df['filter_passed_pulse'] = df['filter_passed_pulse'].astype(bool)
        return  df[(df['filter_passed_L']) | (df['filter_passed_pulse'])]
    
    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on precursor quantity')
        df = df.copy()
        # cols = ['precursor_quantity_L', 'precursor_translated_L', 'ms1_translated_L', 'precursor_translated_pulse', 'ms1_translated_pulse', 'precursor_quantity_L', 'precursor_quantity_pulse']
        
        df.loc[:, 'Precursor.Quantity'] = df['precursor_quantity_L'].fillna(0) + df['precursor_quantity_pulse'].fillna(0)
        
        df.loc[:, 'precursor_quantity_pulse_L_ratio'] = df['precursor_quantity_pulse'] / df['precursor_quantity_L'] 
       
        
        df.loc[:, 'precursor_quantity_pulse_L_ratio'] = df['precursor_quantity_pulse_L_ratio'].replace([np.inf, 0], np.nan)
        
        
        df.loc[:, 'Lib.PG.Q.Value'] = 0
        
        return df
    
    def compute_protein_level_ratios(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
    
        df = df.dropna(subset=['precursor_quantity_pulse_L_ratio'])
    
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
            
            def combined_median(pre_quantity):
                
                if len(pre_quantity.dropna()) <= 1:  # Remove NaNs before counting
                    return np.nan
                else:
                    
                    pre_quantity = np.log2(pre_quantity)  # Log-transform the combined series
                    return 2**np.median(pre_quantity)  # Return the median of the log-transformed values
    
            # Group by protein group and apply the custom aggregation
            grouped_run = run_df.groupby(['protein_group']).apply(lambda x: pd.Series({
                'pulse_L_ratio': combined_median(x['precursor_quantity_pulse_L_ratio'])
            })).reset_index()
    
            grouped_run['Run'] = run
            runs_list.append(grouped_run)
    
        result = pd.concat(runs_list, ignore_index=True)
        result['channel'] = 'ratio'
        return result
    
    def get_protein_level_channel_mask(self, df):
        grouped_dfs = dict(tuple(df.groupby(['Run', 'protein_group'])))
        count_passes = 0
        # Print each group
        protein_groups_channel_list = []
        channel_data = []
        channel_df = pd.DataFrame()
        
        for key, sub_df in grouped_dfs.items():
            
            if np.sum(np.isfinite(sub_df['precursor_quantity_pulse_L_ratio'])) >= 2:
               
                count_passes += 1
            else:
                is_light = False
                is_pulse = False
                is_valid = False
                
                if np.sum(np.isfinite(sub_df['precursor_quantity_L'])) >= 2:
                    is_light = True
                    is_valid = True

        
                if np.sum(np.isfinite(sub_df['precursor_quantity_pulse'])) >= 2:
                    is_pulse = True
                    is_valid = True
                    
                if is_pulse and is_light:
                    is_valid = False

                if is_valid:
                    if is_light:
                        #print(f"\nGroup {key}: is light")
                        channel_data = {'Run':[key[0]],
                                       'protein_group': [key[1]],
                                       'channel':['L']}
                        
                        channel_df = pd.DataFrame(channel_data)
                        
                    elif is_pulse:
                        #print(f"\nGroup {key}: is pulse")
                        channel_data = {'Run':[key[0]],
                                       'protein_group': [key[1]],
                                       'channel':['pulse']}
                        
                        channel_df = pd.DataFrame(channel_data)
                    
        
                    protein_groups_channel_list.append(channel_df)
        
        
        channel_df = pd.concat(protein_groups_channel_list, axis=0)
        channel_df = channel_df.groupby(['protein_group', 'Run'], as_index=False).first()
        channel_df
        return channel_df

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
    
    def merge_data(self, protein_group_ratios, protein_channel_mask, protein_intensities_dlfq):
        # protein_groups = pd.merge(protein_group_ratios, protein_intensites_unnormalized,on=['protein_group','Run'], how='left')
        protein_groups = pd.merge(protein_group_ratios, protein_channel_mask,on=['protein_group','Run', 'channel'], how='outer')

        protein_groups = pd.merge(protein_groups, protein_intensities_dlfq,on=['protein_group','Run'], how='left')
        protein_groups = protein_groups.dropna(subset=['normalized_intensity'])
        
        protein_groups['L'] = 0
        protein_groups['pulse'] = 0
        print(protein_groups.columns.values.tolist())
        return protein_groups
    
    def extract_M_and_L(self, df):
        df.loc[df['channel'] == 'ratio', 'L'] = df['normalized_intensity'] / (df['pulse_L_ratio'] + 1)
        df.loc[df['channel'] == 'ratio', 'pulse'] = df['normalized_intensity'] - df['L']
        df.loc[df['channel'] == 'L', 'L'] = df['normalized_intensity']
        df.loc[df['channel'] == 'pulse', 'pulse'] = df['normalized_intensity']
        df.replace(0, np.nan, inplace=True)
        
        return df
    

    