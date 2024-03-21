# -*- coding: utf-8 -*-
"""
Created on Thu Jan 25 10:49:28 2024

@author: robbi
"""
import pandas as pd
import numpy as np
import time
from tqdm import tqdm
from icecream import ic
from .utils import manage_directories

from silac_dia_tools.pipeline.utils import dlfq_functions as dlfq


class DynamicDiaSis:
    def __init__(self, path, filtered_report):
        self.path = path
        self.filtered_report = filtered_report
        self.update = True
        
        self.formatted_precursors = None
        self.protein_groups = None
        
    def generate_protein_groups(self):
        start_time = time.time()
        # formatting SILAC channels
        self.formatted_precursors = self.format_silac_channels(self.filtered_report)
        
        # Calculate precursor ratios
        self.formatted_precursors = self.calculate_precursor_ratios(self.formatted_precursors)
        
        # Calculate global heavy reference df 
        self.href_df = self.calculate_href_intensities(self.formatted_precursors)
        
        # Calculate protein level ratios
        self.protein_groups = self.compute_protein_level(self.formatted_precursors)
        
        # Merge href with protein groups to generate normalized intensities         
        self.protein_groups = self.href_normalization(self.protein_groups, self.href_df)
        
        # Format and output protein groups for light, NSP, and href dfs
        self.output_protein_groups(self.protein_groups, self.path)
        
        end_time = time.time()
        print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        return self.formatted_precursors, self.protein_groups
    
    def format_silac_channels(self, df):
        print('Formatting SILAC channels')
        # Pivot for each label
        pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' L')
        pivot_M = df[df['Label'] == 'M'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' M')
        pivot_H = df[df['Label'] == 'H'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' H')
        
        # Merge the pivoted DataFrames
        formatted_precursors = pd.concat([pivot_L, pivot_M, pivot_H], axis=1)
        
        # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        formatted_precursors.reset_index(inplace=True)
        
        # Replace 0, inf, -inf with NaN for the specified columns
        formatted_precursors['Precursor.Translated H'] = formatted_precursors['Precursor.Translated H'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Precursor.Translated L'] = formatted_precursors['Precursor.Translated L'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Precursor.Translated M'] = formatted_precursors['Precursor.Translated M'].replace([0, np.inf, -np.inf], np.nan)
        
        formatted_precursors['Ms1.Translated H'] = formatted_precursors['Ms1.Translated H'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Ms1.Translated L'] = formatted_precursors['Ms1.Translated L'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Ms1.Translated L'] = formatted_precursors['Ms1.Translated M'].replace([0, np.inf, -np.inf], np.nan)
 
        return formatted_precursors    
    
    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        # Calculate ratios for all chanels (Precursor.Quantity is the total intensity of all 3 chanels, the default diann value has been overwritten at this point)
        
        df['Precursor.Translated L/H'] = df['Precursor.Translated L'] / df['Precursor.Translated H']
        df['Ms1.Translated L/H'] = df['Ms1.Translated L'] / df['Ms1.Translated H']
        
        df['Precursor.Translated M/H'] = df['Precursor.Translated M'] / df['Precursor.Translated H'] 
        df['Ms1.Translated M/H'] = df['Ms1.Translated M'] / df['Ms1.Translated H']
        
        return df
    
    def calculate_href_intensities(self, df):
        df = df.copy(deep = True)
        df = df.dropna(subset=['Precursor.Translated H','Ms1.Translated H'])
        
        def combined_median(ms1_series, precursor_series):
            combined_series = np.concatenate([ms1_series, precursor_series])
            combined_series = np.log10(combined_series)  # Log-transform the combined series
            return np.median(combined_series)  # Return the median of the log-transformed values
         
        # Group by protein group and apply the custom aggregation
        grouped = df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
            'href': combined_median(x['Ms1.Translated H'], x['Precursor.Translated H']) 
        })).reset_index()
       
        return grouped[['Protein.Group', 'href']]
    
    def compute_protein_level(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
        
        # Drop non shared precursor and ms1 translated ratios for L/H and M/H
        df = df.dropna(subset=['Precursor.Translated L/H','Ms1.Translated L/H'])
        df = df.dropna(subset=['Precursor.Translated M/H','Ms1.Translated M/H'])
    
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
    
            def combined_median(ms1_series, precursor_series):
                combined_series = np.concatenate([ms1_series, precursor_series])
                combined_series = np.log10(combined_series)  # Log-transform the combined series
                return np.median(combined_series)  # Return the median of the log-transformed values
             
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
                'L/H ratio': combined_median(x['Ms1.Translated L/H'], x['Precursor.Translated L/H']),
                'M/H ratio': combined_median(x['Ms1.Translated M/H'], x['Precursor.Translated M/H'])})).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
    
        result = pd.concat(runs_list, ignore_index=True)

        cols = ['Run','Protein.Group', 'L/H ratio', 'M/H ratio']
    
        # Returning the dataframe with specified columns
        return result[cols]   
    
    # Adjust unnormalized intensities
    def href_normalization(self, protein_groups, href):
        print('Calculating adjusted intensities using reference')
        # Merge the href_df onto protein groups containing optimized ratios
        merged_df = protein_groups.merge(href, on='Protein.Group', how='left')
        
        # Obtain normalized light intensities by adding the L/H ratio to the heavy refference in log space
        merged_df['L_norm'] = merged_df['L/H ratio'] + merged_df['href']
        merged_df['M_norm'] = merged_df['M/H ratio'] + merged_df['href']
        
        # reverse log data to output protein intensities
        merged_df['L_norm'] = 10**merged_df['L_norm'] 
        merged_df['M_norm'] = 10**merged_df['M_norm']
        merged_df['href'] = 10**merged_df['href']
        return merged_df
    
    def output_protein_groups(self, df, path):
        manage_directories.create_directory(self.path, 'protein_groups')
        print(f'Outputing normalized protein intensities to {path}/protein_groups')
        cols = ['Run', 'Protein.Group', 'href', 'M_norm', 'L_norm']
        df = df[cols]
        df = df.rename(columns={'href': 'H', 'M_norm': 'M', 'L_norm': 'L'})

        # Pivoting for 'H'
        h_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='H')
        
        # Pivoting for 'M'
        m_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='M')
        
        # Pivoting for 'L'
        l_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='L')
        
        # then output each table to csv for h.href, l.href, m.href
        h_pivot_df.to_csv(f'{path}/protein_groups/href.csv', sep=',')
        m_pivot_df.to_csv(f'{path}/protein_groups/nsp.csv', sep=',')
        l_pivot_df.to_csv(f'{path}/protein_groups/light.csv', sep=',')

        return h_pivot_df, m_pivot_df, l_pivot_df



class DiaSis:
    def __init__(self, path, filtered_report):
        self.path = path
        self.filtered_report = filtered_report
        self.update = True
        
        self.formatted_precursors = None
        self.protein_groups = None
        self.href_df = None
        
    def generate_protein_groups(self):
        start_time = time.time()
        # formatting channels
        self.formatted_precursors = self.format_silac_channels(self.filtered_report)
        
        # calculate L/H ratios for each precursor       
        self.formatted_precursors = self.calculate_precursor_ratios(self.formatted_precursors)
        
        # Calculate global href df
        self.href_df = self.calculate_precursor_href_intensities(self.formatted_precursors)
        
        # Compute protein level ratios
        self.protein_groups = self.compute_protein_level(self.formatted_precursors)
        
        # Merge href df with protein level ratios to normalize light intensities
        self.protein_groups = self.href_normalization(self.protein_groups, self.href_df) 
        
        # Output protein groups table of Light intensities and corrosponding href df
        self.output_protein_groups(self.protein_groups, self.path)
        
        end_time = time.time()
        print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        
        return self.formatted_precursors, self.protein_groups
    
    def format_silac_channels(self, df):
        print('Formatting SILAC channels')
        # Pivot for each label
        pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' L')
        pivot_H = df[df['Label'] == 'H'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' H')
        
        # Merge the pivoted DataFrames
        formatted_precursors = pd.concat([pivot_L, pivot_H], axis=1)
        
        # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        formatted_precursors.reset_index(inplace=True)
        
        # Replace 0, inf, -inf with NaN for the specified columns
        formatted_precursors['Precursor.Translated H'] = formatted_precursors['Precursor.Translated H'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Precursor.Translated L'] = formatted_precursors['Precursor.Translated L'].replace([0, np.inf, -np.inf], np.nan)
        
        formatted_precursors['Ms1.Translated H'] = formatted_precursors['Ms1.Translated H'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Ms1.Translated L'] = formatted_precursors['Ms1.Translated L'].replace([0, np.inf, -np.inf], np.nan)
        
        return formatted_precursors

    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        df['Precursor.Translated L/H'] = df['Precursor.Translated L'] / df['Precursor.Translated H'] 
        df['Ms1.Translated L/H'] = df['Ms1.Translated L'] / df['Ms1.Translated H'] 
      
        return df

    def calculate_precursor_href_intensities(self, df):
        print('Calculate href df')
        df = df.copy(deep = True)
        df = df.dropna(subset=['Precursor.Translated H','Ms1.Translated H'])
        
        def combined_median(ms1_series, precursor_series):
            combined_series = np.concatenate([ms1_series, precursor_series])
            combined_series = np.log10(combined_series)  # Log-transform the combined series
            return np.median(combined_series)  # Return the median of the log-transformed values
       
        # Group by protein group and apply the custom aggregation
        grouped = df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
            'href': combined_median(x['Ms1.Translated H'], x['Precursor.Translated H']) 
        })).reset_index()
       
        return grouped[['Protein.Group', 'href']]
    
   
    def compute_protein_level(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
        
        # Drop non shared precursor and ms1 translated ratios
        df = df.dropna(subset=['Precursor.Translated L/H','Ms1.Translated L/H'])
        
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
    
            def combined_median_ratios(ms1_series, precursor_series):
                combined_series = np.concatenate([ms1_series, precursor_series])
                combined_series = np.log10(combined_series)  # Log-transform the combined series
                return np.median(combined_series)  # Return the median of the log-transformed values
    
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby('Protein.Group').apply(lambda x: pd.Series({
                'L/H ratio': combined_median_ratios(x['Ms1.Translated L/H'], x['Precursor.Translated L/H'])
            })).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
    
        result = pd.concat(runs_list, ignore_index=True)
        cols = ['Run', 'Protein.Group', 'L/H ratio']
        return result[cols]

    
    def href_normalization(self, protein_groups, href):
        print('Calculating adjusted intensities using reference')
        # Merge the href_df onto protein groups containing optimized ratios
        merged_df = protein_groups.merge(href, on='Protein.Group', how='inner')
        
        # Obtain normalized light intensities by adding the L/H ratio to the heavy refference in log space
        merged_df['L_norm'] = merged_df['L/H ratio'] + merged_df['href']
        
        return merged_df
    
    
    def output_protein_groups(self, df, path):
        manage_directories.create_directory(self.path, 'protein_groups')
        print(f'Outputing normalized protein intensities to {path}/protein_groups')
        
        # Subset and rename columns
        df = df[['Run', 'Protein.Group', 'href', 'L_norm']]
        df = df.rename(columns={'href': 'H', 'L_norm': 'L'})

        # Pivoting for 'H' to produce href output in wide format
        h_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='H')
        
        # Pivoting for 'L' to produce normalized light intensites in wide format
        l_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='L')

        # then output each table to csv 
        h_pivot_df.to_csv(f'{path}/protein_groups/href.csv', sep=',')
        l_pivot_df.to_csv(f'{path}/protein_groups/light.csv', sep=',')

        return h_pivot_df, l_pivot_df


class DynamicSilac:    
    def __init__(self, path, filtered_report, pulse_channel):
        self.path = path
        self.filtered_report = filtered_report
        self.update = True
        self.pulse_channel = pulse_channel
        self.formatted_precursors = None
        self.protein_groups = None
        
    def generate_protein_groups(self):  
        start_time = time.time()        
        # formatting and ratios
        self.formatted_precursors = self.format_silac_channels(self.filtered_report)
        self.formatted_precursors = self.calculate_precursor_ratios(self.formatted_precursors)
        self.protein_groups = self.compute_protein_level_ratios(self.formatted_precursors)
        
        # Adjusting intensities and outputing data
        dlfq = self.perform_lfq(self.path)
        self.protein_groups = self.merge_dlfq_intensities(self.protein_groups, dlfq)
        self.output_protein_groups(self.protein_groups, self.path)
        end_time = time.time()
        print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        return self.formatted_precursors, self.protein_groups
        
    def format_silac_channels(self, df):
        # Pivot for each label
        pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' L')
        pivot_M = df[df['Label'] == f'{self.pulse_channel}'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(f' {self.pulse_channel}')
        
        # Merge the pivoted DataFrames
        formatted_precursors = pd.concat([pivot_L, pivot_M], axis=1)
        # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        formatted_precursors.reset_index(inplace=True)
    
        # Replace 0, inf, -inf with NaN for the specified columns
        formatted_precursors[f'Precursor.Translated {self.pulse_channel}'] = formatted_precursors[f'Precursor.Translated {self.pulse_channel}'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Precursor.Translated L'] = formatted_precursors['Precursor.Translated L'].replace([0, np.inf, -np.inf], np.nan)
        
        formatted_precursors[f'Ms1.Translated {self.pulse_channel}'] = formatted_precursors[f'Ms1.Translated {self.pulse_channel}'].replace([0, np.inf, -np.inf], np.nan)
        formatted_precursors['Ms1.Translated L'] = formatted_precursors['Ms1.Translated L'].replace([0, np.inf, -np.inf], np.nan)
        return formatted_precursors
    
    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        # Calculate ratios for all chanels (Precursor.Quantity is the total intensity of all 3 chanels, the default diann value has been overwritten at this point)
        df['Precursor.Quantity'] = df['Precursor.Translated L'].fillna(0) +  df[f'Precursor.Translated {self.pulse_channel}'].fillna(0)
        
        df[f'Precursor.Translated {self.pulse_channel}/L'] = df[f'Precursor.Translated {self.pulse_channel}'] / df['Precursor.Translated L']
        df[f'Ms1.Translated {self.pulse_channel}/L'] = df[f'Ms1.Translated {self.pulse_channel}'] / df['Ms1.Translated L']
        # Replace invalid values with NaN and drop them
        df[f'Precursor.Translated {self.pulse_channel}/L'] = df[f'Precursor.Translated {self.pulse_channel}/L'].replace([0, np.inf, -np.inf], np.nan)
        df[f'Ms1.Translated {self.pulse_channel}/L'] = df[f'Ms1.Translated {self.pulse_channel}/L'].replace([0, np.inf, -np.inf], np.nan)
        df['Lib.PG.Q.Value'] = 0
        return df
    
    def compute_protein_level_ratios(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
        
        
        df = df.dropna(subset=[f'Ms1.Translated {self.pulse_channel}/L', f'Precursor.Translated {self.pulse_channel}/L'])
        
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
    
            def combined_median(ms1_series, precursor_series):
                combined_series = np.concatenate([ms1_series, precursor_series])
                combined_series = np.log2(combined_series)  # Log-transform the combined series
                return 2 ** np.median(combined_series)  # Return the median of the log-transformed values
           
            def valid_median(series):
                # valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
                return series.median()
    
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
                f'{self.pulse_channel}/L ratio': combined_median(x[f'Ms1.Translated {self.pulse_channel}/L'], x[f'Precursor.Translated {self.pulse_channel}/L'])
            })).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
    
        result = pd.concat(runs_list, ignore_index=True)
                
        cols = ['Run','Protein.Group', f'{self.pulse_channel}/L ratio']
    
        # Returning the dataframe with specified columns
        return result[cols]
    
    def perform_lfq(self, path):
        dlfq_formatted_precursors = self.formatted_precursors
        # ic(self.formatted_precursors)
        manage_directories.create_directory(self.path, 'dlfq')

        dlfq_formatted_precursors.to_csv(f'{path}/dlfq/dflq_formatted_report.tsv', sep='\t')
        dlfq_output_file = f'{path}/dlfq/dlfq_protein_intensities.tsv'
        
        dlfq.run_lfq(f'{path}/dlfq/dflq_formatted_report.tsv', file=dlfq_output_file, num_cores=1)
        dlfq_df = pd.read_csv(dlfq_output_file, sep='\t')
       
        # Drop the 'Unnamed: 0' column
        dlfq_df = dlfq_df.drop(columns=['Unnamed: 0', 'protein'])
       
        # Melt the DataFrame
        long_df = pd.melt(dlfq_df, id_vars=['Protein.Group'], var_name='Run', value_name='Intensity')
        # merged_dlfq = merge_dlfq_intensities(protein_groups_unnormalized, long_df)
        
        return long_df
    
    def merge_dlfq_intensities(self, df, dlfq):     
        # Merge the original DataFrame with the h_ref DataFrame
        merged_df = df.merge(dlfq, on=['Protein.Group','Run'], how='inner')
        
        # Use ratios to calculate intensiteis
        merged_df[f'{self.pulse_channel}_norm'] =  merged_df['Intensity'] / (merged_df[f'{self.pulse_channel}/L ratio' ] + 1)
        merged_df['L_norm'] = merged_df[f'{self.pulse_channel}/L ratio'] *  merged_df[f'{self.pulse_channel}_norm']
     
        return merged_df
    
    def output_protein_groups(self, df, path):
        manage_directories.create_directory(self.path, 'protein_groups')

        cols = ['Run', 'Protein.Group', f'{self.pulse_channel}_norm', 'L_norm']
        df = df[cols]
        df = df.rename(columns={ f'{self.pulse_channel}_norm': f'{self.pulse_channel}', 'L_norm': 'L'})
        
        # Pivoting for 'M'
        m_pivot_df = df.pivot(index='Protein.Group', columns='Run', values=f'{self.pulse_channel}')
        
        # Pivoting for 'L'
        l_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='L')
        
        # then output each table to csv for h.href, l.href, m.href
        m_pivot_df.to_csv(f'{path}protein_groups/nsp.csv', sep=',')
        l_pivot_df.to_csv(f'{path}protein_groups/light.csv', sep=',')
    
        return  m_pivot_df, l_pivot_df
    

      
      
      
      
      