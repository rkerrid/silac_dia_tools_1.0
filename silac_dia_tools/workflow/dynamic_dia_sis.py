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
        
        protein_ratios = pd.merge(L_pg, M_pg, on=['Run', 'protein_group'], how='outer' )
        protein_ratios = pd.merge(href, protein_ratios, on='protein_group', how='outer')
        
        protein_ratios = protein_ratios.rename(columns={'precursor_quantity_H' : 'href'})
        protein_ratios['L'] = protein_ratios['L/H ratio'] * protein_ratios['href']
        protein_ratios['M'] = protein_ratios['M/H ratio'] * protein_ratios['href']
        
        return protein_ratios
        
    def get_filtered_set(self, df):
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
        
        return href_df
    
    def get_precursor_ratios(self, precursors):
        precursors['precursor_translated_LH'] = precursors['precursor_translated_L'] / precursors['precursor_translated_H'] 
        precursors['ms1_translated_LH'] = precursors['ms1_translated_L'] / precursors['ms1_translated_H'] 
        precursors['precursor_translated_MH'] = precursors['precursor_translated_M'] / precursors['precursor_translated_H'] 
        precursors['ms1_translated_MH'] = precursors['ms1_translated_M'] / precursors['ms1_translated_H']
        
        columns_to_replace = ['precursor_translated_LH', 'ms1_translated_LH', 'precursor_translated_MH', 'ms1_translated_MH']
        precursors[columns_to_replace] = precursors[columns_to_replace].replace([np.inf, -np.inf, 0.0], np.nan)
        
        index_cols = ['Run','protein_group','precursor_id']
        LH_precursors = precursors[index_cols + ['precursor_translated_LH', 'ms1_translated_LH']]
        LH_precursors = LH_precursors.dropna(subset=['precursor_translated_LH', 'ms1_translated_LH'])
        MH_precursors = precursors[index_cols + ['precursor_translated_MH', 'ms1_translated_MH']]
        MH_precursors = MH_precursors.dropna(subset=['precursor_translated_MH', 'ms1_translated_MH'])
        
        return LH_precursors, MH_precursors 

    def get_protein_ratios(self, df, channel):
        runs = df['Run'].unique()
        runs_list = []
        
        for run in runs:
            run_df = df[df['Run'] == run]
        
            def combined_median(ms1_series, precursor_series):
                combined_series = np.concatenate([ms1_series, precursor_series])
                combined_series = np.log2(combined_series)  # Log-transform the combined series
                return 2**np.median(combined_series)  # Return the median of the log-transformed values
             
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby(['protein_group']).apply(lambda x: pd.Series({
                f'{channel}/H ratio': combined_median(x[f'ms1_translated_{channel}H'], x[f'precursor_translated_{channel}H'])})).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
        
        result = pd.concat(runs_list, ignore_index=True)
        
        return result

# class DynamicDiaSis_:
#     def __init__(self, path, filtered_report):
#         self.path = path
#         self.filtered_report = filtered_report
#         self.update = True
        
#         self.formatted_precursors = None
#         self.protein_groups = None
        
#     def generate_protein_groups(self):
#         start_time = time.time()
        
#         # formatting SILAC channels
#         h_precursors_df, LH_df, MH_df  = self.format_silac_channels(self.filtered_report)
      
#         # Calculate global heavy reference df 
#         href_df = self.calculate_href_intensities(h_precursors_df)
        
#         # calculate protein level ratios
#         LH_protein_df = self.compute_protein_level_ratios(LH_df, 'L')
#         MH_protein_df = self.compute_protein_level_ratios(MH_df, 'M')
        
#         # Merge href with protein groups to generate normalized intensities         
#         LH_protein_df, MH_protein_df = self.href_normalization( LH_protein_df, MH_protein_df, href_df)
        
#         # oputput data
#         self.output_protein_groups( LH_protein_df, MH_protein_df, href_df, self.path)
        

#         end_time = time.time()
#         print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        
    
#     def format_silac_channels(self, df):
#         # pivot table
#         df = df.pivot_table(index=['Run','Protein.Group', 'Precursor.Id'], columns='Label', values = ['Ms1.Translated', 'Precursor.Quantity', 'Precursor.Translated'])
        
#         # format href, LH and MH precursors
#         href_df = self.format_reference_silac_channels(df)
#         LH_df = self.format_LH_silac_channels(df)
#         MH_df = self.format_MH_silac_channels(df)
        
#         # calculate ratios for LH and MH dfs
#         LH_df = self.calculate_precursor_ratios(LH_df, 'L')
#         MH_df = self.calculate_precursor_ratios(MH_df, 'M')

#         return href_df, LH_df, MH_df    
    
#     def format_reference_silac_channels(self, df):
#         # Access the 'Ms1.Translated' and 'Precursor.Translated' columns under the 'H' label
#         ms1_translated_H = df.loc[:, ('Ms1.Translated', 'H')]
#         precursor_translated_H = df.loc[:, ('Precursor.Translated', 'H')]
        
#         # Combine into a new DataFrame
#         combined_df = pd.DataFrame({
#             'Ms1.Translated_H': ms1_translated_H,
#             'Precursor.Translated_H': precursor_translated_H
#         })
        
#         # Reset index to include 'Run', 'Protein.Group', and 'Precursor.Id' as columns
#         combined_H_df = combined_df.reset_index()
        
#         # Rename columns for clarity
#         combined_H_df.columns = ['Run', 'Protein.Group', 'Precursor.Id', 'Ms1.Translated_H', 'Precursor.Translated_H']
        
#         # drop rows where there are no valid pairs of Precursor and Ms1 translated
#         combined_H_df = combined_H_df.dropna(subset=['Precursor.Translated_H','Ms1.Translated_H'])
        
#         return combined_H_df
        
#     def format_LH_silac_channels(self, df):
#         # Access the 'Ms1.Translated' and 'Precursor.Translated' columns under the 'H' and 'L' label
#         ms1_translated_L = df.loc[:, ('Ms1.Translated', 'L')]
#         precursor_translated_L = df.loc[:, ('Precursor.Translated', 'L')]
#         ms1_translated_H = df.loc[:, ('Ms1.Translated', 'H')]
#         precursor_translated_H = df.loc[:, ('Precursor.Translated', 'H')]
        
#         # Combine into a new DataFrame
#         combined_df = pd.DataFrame({
#             'Ms1.Translated_H': ms1_translated_H,
#             'Precursor.Translated_H': precursor_translated_H,
#             'Ms1.Translated_L': ms1_translated_L,
#             'Precursor.Translated_L': precursor_translated_L
#         })
        
#         # Reset index to include 'Run', 'Protein.Group', and 'Precursor.Id' as columns
#         combined_LH_df = combined_df.reset_index()
        
#         # Rename columns for clarity
#         combined_LH_df.columns = ['Run', 'Protein.Group', 'Precursor.Id', 'Ms1.Translated_H', 'Precursor.Translated_H', 'Ms1.Translated_L', 'Precursor.Translated_L']
               
#         return combined_LH_df
    
#     def format_MH_silac_channels(self, df):
#         # Access the 'Ms1.Translated' and 'Precursor.Translated' columns under the 'H' and 'M' label
#         ms1_translated_M = df.loc[:, ('Ms1.Translated', 'M')]
#         precursor_translated_M = df.loc[:, ('Precursor.Translated', 'M')]
#         ms1_translated_H = df.loc[:, ('Ms1.Translated', 'H')]
#         precursor_translated_H = df.loc[:, ('Precursor.Translated', 'H')]
        
#         # Combine into a new DataFrame
#         combined_df = pd.DataFrame({
#             'Ms1.Translated_H': ms1_translated_H,
#             'Precursor.Translated_H': precursor_translated_H,
#             'Ms1.Translated_M': ms1_translated_M,
#             'Precursor.Translated_M': precursor_translated_M
#         })
        
#         # Reset index to include 'Run', 'Protein.Group', and 'Precursor.Id' as columns
#         combined_MH_df = combined_df.reset_index()
        
#         # Rename columns for clarity
#         combined_MH_df.columns = ['Run', 'Protein.Group', 'Precursor.Id', 'Ms1.Translated_H', 'Precursor.Translated_H', 'Ms1.Translated_M', 'Precursor.Translated_M']
  
#         return combined_MH_df
        
#     def calculate_precursor_ratios(self, df, channel):
#         print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
#         # Calculate ratios for channel        
#         df[f'Precursor.Translated {channel}/H'] = df[f'Precursor.Translated_{channel}'] / df['Precursor.Translated_H']
#         df[f'Ms1.Translated {channel}/H'] = df[f'Ms1.Translated_{channel}'] / df['Ms1.Translated_H']
        
#         columns_to_replace = [f'Precursor.Translated {channel}/H', f'Ms1.Translated {channel}/H']
#         df[columns_to_replace] = df[columns_to_replace].replace([np.inf, -np.inf, 0.0], np.nan)
        
#         df = df.dropna(subset=[f'Precursor.Translated {channel}/H',f'Ms1.Translated {channel}/H'])
        
#         return df
  
#     def calculate_href_intensities(self, df):
        
#         def combined_median(ms1_series, precursor_series):
#             combined_series = np.concatenate([ms1_series, precursor_series])
#             combined_series = np.log10(combined_series)  # Log-transform the combined series
#             return np.median(combined_series)  # Return the median of the log-transformed values
         
#         # Group by protein group and apply the custom aggregation
#         grouped = df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
#             'href': combined_median(x['Ms1.Translated_H'], x['Precursor.Translated_H']) 
#         })).reset_index()
        
        
#         return  grouped[['Protein.Group', 'href']]
    
#     def compute_protein_level_ratios(self, df, channel):
#         runs = df['Run'].unique()
#         runs_list = []
        
#         for run in tqdm(runs, desc=f'Computing protein level ratios for each run, channel: {channel}'):
#             run_df = df[df['Run'] == run]
        
#             def combined_median(ms1_series, precursor_series):
#                 combined_series = np.concatenate([ms1_series, precursor_series])
#                 # print(combined_series)
#                 combined_series = np.log10(combined_series)  # Log-transform the combined series
#                 return np.median(combined_series)  # Return the median of the log-transformed values
             
#             # Group by protein group and apply the custom aggregation
#             grouped = run_df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
#                 f'{channel}/H ratio': combined_median(x[f'Ms1.Translated {channel}/H'], x[f'Precursor.Translated {channel}/H'])})).reset_index()
            
#             grouped['Run'] = run
#             runs_list.append(grouped)
        
#         result = pd.concat(runs_list, ignore_index=True)
        
#         cols = ['Run','Protein.Group', f'{channel}/H ratio']
         
#         return result[cols] 
    
#     def href_normalization(self, LH_protein_df, MH_protein_df, href_df):
#         # Merge the href_df onto protein groups containing optimized ratios
#         merged_df_LH = LH_protein_df.merge(href_df, on='Protein.Group', how='left')
#         merged_df_MH = MH_protein_df.merge(href_df, on='Protein.Group', how='left')
        
#         # Obtain normalized light intensities by adding the L/H ratio to the heavy refference in log space
#         merged_df_LH['L_norm'] = merged_df_LH['L/H ratio'] + merged_df_LH['href']
#         merged_df_MH['M_norm'] = merged_df_MH['M/H ratio'] + merged_df_MH['href']
        
#         # reverse log data to output protein intensities
#         merged_df_LH['L_norm'] = 10**merged_df_LH['L_norm'] 
#         merged_df_MH['M_norm'] = 10**merged_df_MH['M_norm']
#         href_df['href'] = 10**href_df['href']
        
#         return merged_df_LH, merged_df_MH    
    
#     def output_protein_groups(self, LH_protein_df, MH_protein_df, href_df, path):
#         manage_directories.create_directory(self.path, 'protein_groups')
#         print(f'Outputing normalized protein intensities to {path}/protein_groups')
#         LH_protein_df = LH_protein_df.rename(columns={'L_norm': 'L'})
#         MH_protein_df = MH_protein_df.rename(columns={'M_norm': 'M'})
        
#         # Pivoting for 'L'
#         l_pivot_df = LH_protein_df.pivot(index='Protein.Group', columns='Run', values='L')
        
#         # Pivoting for 'M'
#         m_pivot_df = MH_protein_df.pivot(index='Protein.Group', columns='Run', values='M')

#         # then output each table to csv for h.href, l.href, m.href
#         # h_pivot_df.to_csv(f'{path}/protein_groups/href.csv', sep=',')
#         href_df.to_csv(f'{path}/protein_groups/href.csv', sep=',')
#         m_pivot_df.to_csv(f'{path}/protein_groups/nsp.csv', sep=',')
#         l_pivot_df.to_csv(f'{path}/protein_groups/light.csv', sep=',')