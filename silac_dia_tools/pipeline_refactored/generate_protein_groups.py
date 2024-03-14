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


class DynamicSilacDiaSis:
    def __init__(self, path, filtered_report):
        self.path = path
        self.filtered_report = filtered_report
        self.update = True
        
        self.formatted_precursors = None
        self.protein_groups = None
        
    def generate_protein_groups(self):
        start_time = time.time()
        # formatting and ratios
        self.formatted_precursors = self.format_silac_channels(self.filtered_report)
        self.formatted_precursors = self.calculate_precursor_ratios(self.formatted_precursors)
        self.protein_groups = self.compute_protein_level(self.formatted_precursors)
        # Adjusting intensities and outputing data
        self.href_df = self.calculate_href_intensities(self.formatted_precursors)
        self.protein_groups = self.href_normalization(self.protein_groups, self.href_df)
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
        merged_df = pd.concat([pivot_L, pivot_M, pivot_H], axis=1)
        # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        merged_df.reset_index(inplace=True)
    
        # remove all rows that contain a NaN under the Label H column (i.e., no H precursor is present for that row)
        # Apply dropna on merged_df instead of df
        merged_df = merged_df.dropna(subset=['Precursor.Translated H'])
        # replace precursor quantity with summed silac channels as input for direct lefq and as 'total intensity' for href quantification
        merged_df['Precursor.Quantity'] = merged_df['Precursor.Translated H'] + merged_df['Precursor.Translated M'] + merged_df['Precursor.Translated L'] 
 
        return merged_df
    
    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        # Calculate ratios for all chanels (Precursor.Quantity is the total intensity of all 3 chanels, the default diann value has been overwritten at this point)
        
        df['Precursor.Translated L/H'] = df['Precursor.Translated L'] / df['Precursor.Translated H']
        df['Ms1.Translated L/H'] = df['Ms1.Translated L'] / df['Ms1.Translated H']
        
        df['Precursor.Translated M/H'] = df['Precursor.Translated M'] / df['Precursor.Translated H'] 
        df['Ms1.Translated M/H'] = df['Ms1.Translated M'] / df['Ms1.Translated H']
        
        return df
    
    def calculate_href_intensities(self, df):
        
        def combined_median(ms1_series, precursor_series):
            # Replace invalid values with NaN and drop them
            valid_ms1 = ms1_series.replace([0, np.inf, -np.inf], np.nan).dropna()
            valid_precursor = precursor_series.replace([0, np.inf, -np.inf], np.nan).dropna()
       
            # Ensure at least 3 valid values in each series before combining
            if len(valid_ms1) >= 1 or len(valid_precursor) >= 1:
                combined_series = np.concatenate([valid_ms1, valid_precursor])
                combined_series = np.log10(combined_series)  # Log-transform the combined series
                return np.median(combined_series)  # Return the median of the log-transformed values
            else:
                return np.nan
       
        # Group by protein group and apply the custom aggregation
        grouped = df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
            'href': combined_median(x['Ms1.Translated H'], x['Precursor.Translated H']) 
        })).reset_index()
       
        return grouped[['Protein.Group', 'href']]
    
    def compute_protein_level(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
    
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
    
            def combined_median(ms1_series, precursor_series):
                # Replace invalid values with NaN and drop them
                valid_ms1 = ms1_series.replace([0, np.inf, -np.inf], np.nan).dropna()
                valid_precursor = precursor_series.replace([0, np.inf, -np.inf], np.nan).dropna()
    
                # Ensure at least 3 valid values in each series before combining
                if len(valid_ms1) >= 1 or len(valid_precursor) >= 1:
                    combined_series = np.concatenate([valid_ms1, valid_precursor])
                    combined_series = np.log10(combined_series)  # Log-transform the combined series
                    return np.median(combined_series)  # Return the median of the log-transformed values
                else:
                    return np.nan
    
    
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
        
    def generate_protein_groups(self):
        start_time = time.time()

        # formatting and ratios
        self.formatted_precursors = self.format_silac_channels(self.filtered_report)
        self.formatted_precursors = self.calculate_precursor_ratios(self.formatted_precursors)
        self.protein_groups = self.compute_protein_level(self.formatted_precursors)
        # self.protein_groups = self.compute_protein_level_test(self.formatted_precursors)

        # Adjusting intensities and outputing data
        self.protein_groups = self.calculate_href_intensities(self.protein_groups)
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
        merged_df = pd.concat([pivot_L, pivot_H], axis=1)
        # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        merged_df.reset_index(inplace=True)
    
        # remove all rows that contain a NaN under the Label H column (i.e., no H precursor is present for that row)
        # Apply dropna on merged_df instead of df
        # merged_df = merged_df.dropna(subset=['Precursor.Translated H'])
        # replace precursor quantity with summed silac channels as input for direct lefq and as 'total intensity' for href quantification
        merged_df['Precursor.Quantity'] = merged_df['Precursor.Translated H'] + merged_df['Precursor.Translated L'] 
 
        return merged_df
    
    def calculate_precursor_ratios(self, df):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        # Calculate ratios for all chanels (Precursor.Quantity is the total intensity of all 3 chanels, the default diann value has been overwritten at this point)
        df['Precursor.Translated H/T'] = df['Precursor.Translated H'] / (df['Precursor.Translated L'] + df['Precursor.Translated H'])
        df['Ms1.Translated H/T'] = df['Ms1.Translated H'] / (df['Ms1.Translated L'] + df['Ms1.Translated H'])
        return df    

    def compute_protein_level(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
    
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
    
            def combined_median(ms1_series, precursor_series):
                # Replace invalid values with NaN and drop them
                valid_ms1 = ms1_series.replace([0, np.inf, -np.inf], np.nan).dropna()
                valid_precursor = precursor_series.replace([0, np.inf, -np.inf], np.nan).dropna()
    
                # Ensure at least 3 valid values in each series before combining
                if len(valid_ms1) >= 2 and len(valid_precursor) >= 2:
                    combined_series = np.concatenate([valid_ms1, valid_precursor])
                    combined_series = np.log2(combined_series)  # Log-transform the combined series
                    return 2 ** np.median(combined_series)  # Return the median of the log-transformed values
                else:
                    return np.nan
    
            def valid_sum(series):
                valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
                return valid_series.sum()
    
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
                'H/T ratio': combined_median(x['Ms1.Translated H/T'], x['Precursor.Translated H/T']),
                'Precursor.Quantity': valid_sum(x['Precursor.Quantity'])
            })).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
    
        result = pd.concat(runs_list, ignore_index=True)
        result['H'] = result['H/T ratio']*result['Precursor.Quantity']
        result['L'] = result['Precursor.Quantity'] - result['H']
        cols = ['Run','Protein.Group', 'H', 'L', 'H/T ratio', 'Precursor.Quantity']
    
        # Returning the dataframe with specified columns
        return result[cols]

    # Adjust unnormalized intensities
    def calculate_href_intensities(self, df):
        print('Calculating adjusted intensities using reference')
        df_copy = df.copy(deep=True)
        
        # Calculate median H value and reset index to make it a DataFrame
        h_ref = df_copy.groupby('Protein.Group')['H'].median().reset_index()
        
        # Rename the median column to 'h_ref'
        h_ref = h_ref.rename(columns={'H': 'h_ref'})
        
        # Merge the original DataFrame with the h_ref DataFrame
        merged_df = df.merge(h_ref, on='Protein.Group', how='inner')
        
        # calculate factor to multiply other chanels by dividing href by original H intensity for each PG
        merged_df['factor'] = merged_df['h_ref']/merged_df['H']
        
        # Normalize other chanels with this factor
        merged_df['H_norm'] = merged_df['H']*merged_df['factor']
        merged_df['L_norm'] = merged_df['L']*merged_df['factor']
        
        return merged_df
    
    def output_protein_groups(self, df, path):
        manage_directories.create_directory(self.path, 'protein_groups')
        print(f'Outputing normalized protein intensities to {path}/protein_groups')
        cols = ['Run', 'Protein.Group', 'H_norm', 'L_norm']
        df = df[cols]
        df = df.rename(columns={'H_norm': 'H', 'L_norm': 'L'})

        # Pivoting for 'H'
        h_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='H')
        
        
        # Pivoting for 'L'
        l_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='L')
        
        # then output each table to csv for h.href, l.href, m.href
        h_pivot_df.to_csv(f'{path}/protein_groups/href_href.csv', sep=',')
        l_pivot_df.to_csv(f'{path}/protein_groups/light_href.csv', sep=',')

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
        # Will use Precursor.Translated for quantification
    
        start_time = time.time()
        # Will use Precursor.Translated for quantification
        quantification = 'Precursor.Translated'
        # formatting and ratios
        self.formatted_precursors = self.format_silac_channels(self.filtered_report)
        self.formatted_precursors = self.calculate_precursor_ratios(self.formatted_precursors, quantification)
        self.protein_groups = self.compute_protein_level(self.formatted_precursors)
        # Adjusting intensities and outputing data
        dlfq = self.perform_lfq(self.path)
        self.protein_groups = self.merge_dlfq_intensities(self.protein_groups, dlfq)
        self.output_protein_groups(self.protein_groups, quantification, self.path)
        end_time = time.time()
        print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        return self.formatted_precursors, self.protein_groups
        
    def format_silac_channels(self, df):
        # Pivot for each label
        pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' L')
        pivot_M = df[df['Label'] == f'{self.pulse_channel}'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' M')
        
        # Merge the pivoted DataFrames
        merged_df = pd.concat([pivot_L, pivot_M], axis=1)
        # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        merged_df.reset_index(inplace=True)
    
  
        print("Column 'Label H' not found in DataFrame")
        print(merged_df.columns.values.tolist())
        merged_df['Precursor.Quantity'] = merged_df[f'Precursor.Translated {self.pulse_channel}'] + merged_df['Precursor.Translated L'] 
        return merged_df
    
    def calculate_precursor_ratios(self, df, quantification):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        # Calculate ratios for all chanels (Precursor.Quantity is the total intensity of all 3 chanels, the default diann value has been overwritten at this point)
        df['Precursor.Quantity'] = df['Precursor.Translated L'] +  df['Precursor.Translated {self.pulse_channel}']
        
        df['Precursor.Translated {self.pulse_channel}/T'] = df['Precursor.Translated {self.pulse_channel}'] / df['Precursor.Quantity']
        df['Ms1.Translated {self.pulse_channel}/T'] = df['Ms1.Translated {self.pulse_channel}'] / df['Precursor.Quantity']
        
        df['Precursor.Translated L/T'] = df['Precursor.Translated L'] / df['Precursor.Quantity']
        df['Ms1.Translated L/T'] = df['Ms1.Translated L'] / df['Precursor.Quantity']
        return df
    
    def compute_protein_level(self, df):
        print('Rolling up to protein level')
        runs = df['Run'].unique()
        runs_list = []
    
        for run in tqdm(runs, desc='Computing protein level ratios for each run'):
            run_df = df[df['Run'] == run]
    
            def combined_median(ms1_series, precursor_series):
                # Replace invalid values with NaN and drop them
                valid_ms1 = ms1_series.replace([0, np.inf, -np.inf], np.nan).dropna()
                valid_precursor = precursor_series.replace([0, np.inf, -np.inf], np.nan).dropna()
    
                # Ensure at least 3 valid values in each series before combining
                if len(valid_ms1) >= 2 and len(valid_precursor) >= 2:
                    combined_series = np.concatenate([valid_ms1, valid_precursor])
                    combined_series = np.log2(combined_series)  # Log-transform the combined series
                    return 2 ** np.median(combined_series)  # Return the median of the log-transformed values
                else:
                    return np.nan
    
            def valid_sum(series):
                valid_series = series.replace([0, np.inf, -np.inf], np.nan).dropna()
                return valid_series.sum()
    
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
                '{self.pulse_channel}/T ratio': combined_median(x['Ms1.Translated {self.pulse_channel}/T'], x['Precursor.Translated {self.pulse_channel}/T']),
                'L/T ratio': combined_median(x['Ms1.Translated L/T'], x['Precursor.Translated L/T']),
                'Precursor.Quantity': valid_sum(x['Precursor.Quantity'])
            })).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
    
        result = pd.concat(runs_list, ignore_index=True)
        result['{self.pulse_channel}'] = result['{self.pulse_channel}/T ratio']*result['Precursor.Quantity']
        result['L'] = result['L/T ratio']*result['Precursor.Quantity']
        
        result['Lib.PG.Q.Value'] = 0
        cols = ['Run','Protein.Group', '{self.pulse_channel}', 'L', 'Precursor.Quantity', 'Lib.PG.Q']
    
        # Returning the dataframe with specified columns
        return result[cols]
    
    def perform_lfq(self, path):
        dlfq_formatted_precursors = self.formatted_precursors
        ic(self.formatted_precursors)
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
        merged_df['L_norm'] = merged_df['Precursor.Translated L/T' ] *merged_df['Intensity']
        merged_df[f'{self.pulse_channel}_norm'] = merged_df[f'Precursor.Translated {self.pulse_channel}/T'] *  merged_df['Intensity']
     
        return merged_df
    
    def output_protein_groups(self, df, quantification, path):
        manage_directories.create_directory(self.path, 'protein_groups')

        cols = ['Run', 'Protein.Group', f'{self.pulse_channel}_norm', 'L_norm']
        df = df[cols]
        df = df.rename(columns={ f'{self.pulse_channel}_norm': f'{self.pulse_channel}', 'L_norm': 'L'})
        
        # Pivoting for 'M'
        m_pivot_df = df.pivot(index='Protein.Group', columns='Run', values=f'{self.pulse_channel}')
        
        # Pivoting for 'L'
        l_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='L')
        
        # then output each table to csv for h.href, l.href, m.href
        m_pivot_df.to_csv(f'{path}nsp_dlfq.csv', sep=',')
        l_pivot_df.to_csv(f'{path}light_dlfq.csv', sep=',')
    
        return  m_pivot_df, l_pivot_df
    

      
      
      
      
      