# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 10:14:26 2024

@author: robbi
"""

import pandas as pd
import numpy as np
import operator
import time 

from tqdm import tqdm
import os
from icecream import ic


class Preprocessor:
    def __init__(self, path, pulse_channel,  method, meta_data=None):
        self.path = path
        self.meta_data = meta_data
        self.chunk_size = 180000
     
        self.pulse_channel = pulse_channel
        self.method = method
        
    def import_report(self):
        print('Beginning import report.tsv')
        start_time = time.time()
                
        chunks = []
        contaminants = []
        file_path = f"{self.path}report.tsv"
        count = 1
        # Estimate rows in file size
        file_size_bytes = os.path.getsize(file_path)
        average_row_size_bytes = 1000  # This is an example; you'll need to adjust this based on your data
        # Estimate the number of rows
        estimated_rows = file_size_bytes / average_row_size_bytes
        total_chunks = estimated_rows/self.chunk_size
        for chunk in tqdm(pd.read_table(file_path, sep="\t", chunksize=self.chunk_size), 
                      total=total_chunks, desc='Estimated loading of report.tsv based on file size'):
         
            if self.meta_data is not None:
                chunk = self.subset_based_on_metadata(chunk)
                chunk = self.relabel_run(chunk)
                
            chunk = self.add_label_col(chunk)
            chunk = self.add_passes_filter_col(chunk)
            chunk = self.drop_cols(chunk)
            chunk = self.pivot_data(chunk)
            
            if self.method == 'dynamic_silac_dia':
                chunk = self.format_data_dynamic_silac_dia(chunk)
            elif self.method == 'dynamic_dia_sis':
                chunk = self.format_data_dynamic_dia_sis(chunk)
            elif self.method == 'dia_sis':
                chunk = self.format_data_dia_sis(chunk)
            
            chunk.rename(columns={'Protein.Group':'protein_group','Protein.Ids':'protein_ids', 'Protein.Names':'protein_names', 'Genes':'genes', 'Precursor.Id': 'precursor_id'}, inplace=True)

            chunk, contam_chunk = self.remove_contaminants(chunk)
            
            chunks.append(chunk)
            contaminants.append(contam_chunk)
            
        # append chunks to respective dfs and return  
        filtered_report = pd.concat(chunks, ignore_index=True)
        contaminants_df = pd.concat(contaminants, ignore_index=True)
        print('Finished import')
        end_time = time.time()
        print(f"Time taken for import: {end_time - start_time} seconds")
        
        return filtered_report, contaminants_df


    def subset_based_on_metadata(self, df):       
        filtered_df = df[df['Run'].isin(self.meta_data['Run'])]
        return filtered_df
    
    def relabel_run(self, df):
        run_to_sample = dict(zip(self.meta_data['Run'], self.meta_data['Sample']))
    
        # Apply the mapping to df2['Run'] and raise an error if a 'Run' value doesn't exist in df1
        df['Run'] = df['Run'].map(run_to_sample)
        if df['Run'].isna().any():
            raise ValueError("Some Run values in the report.tsv are not found in the metadata, please ensure metadata is correct.")
            
        return df

    def add_label_col(self, df):
        # Extract the label and add it as a new column
        df['Label'] = df['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]
        
        # Remove the '(SILAC-K|R-([HML]))' part from the 'Precursor.Id' string
        df['Precursor.Id'] = df['Precursor.Id'].str.replace(r'\(SILAC-(K|R)-[HML]\)', '', regex=True)
    
        return df

    def add_passes_filter_col(self, df):
        
        # if data does not pass the following filters set to False
        df['filter_passed'] = (df["Global.PG.Q.Value"] < 0.01) & (df["Precursor.Charge"] > 1) & (df["Channel.Q.Value"] < 0.03)
    
        return df
    
    def drop_cols(self, df):
        # what cols to keep for future workflow
        cols = ['Run',
                 'Protein.Group',
                 'Protein.Ids',
                 'Protein.Names',
                 'Genes',
                 'Precursor.Id',
                 'Precursor.Quantity',
                 'Precursor.Normalised',
                 'Precursor.Translated',
                 'Ms1.Translated',
                 'Label',
                 'filter_passed']
        
        # drop all other cols
        df = df[cols]
        return df


    def remove_contaminants(self, df):
        #chunk_copy = chunk.copy(deep=True)
        contams_mask = df['protein_group'].str.contains('Cont_', case=False, na=False)
        df_filtered = df.loc[~contams_mask].reset_index(drop=True)
        contams = df.loc[contams_mask].reset_index(drop=True)
   
        return df_filtered, contams
    
    
    #pivot data
    def pivot_data(self, df):
        index_cols = ['Run', 'Protein.Group','Protein.Ids','Protein.Names', 'Genes', 'Precursor.Id'] #, 'Protein.Ids', 'Protein.Names', 'Genes', 'Precursor.Id', 'Passes_filter'
        df_p = df.pivot_table(index=index_cols, columns='Label', values = ['Ms1.Translated', 'Precursor.Quantity', 'Precursor.Translated','filter_passed'])
        df_p['filter_passed'] = df_p['filter_passed'].applymap(lambda x: x == 1.0)
        return df_p
    
    def format_data_dynamic_silac_dia(self, df):
        # format new dataframes for future work (3 channels)
        # channel L
        ms1_translated_L = df.loc[:, ('Ms1.Translated', 'L')]
        precursor_translated_L = df.loc[:, ('Precursor.Translated', 'L')]
        precursor_quantity_L = df.loc[:, ('Precursor.Quantity', 'L')]
        passes_filter_L = df.loc[:, ('filter_passed', 'L')]
        
        
        # channel pulse
        ms1_translated_pulse = df.loc[:, ('Ms1.Translated', f'{self.pulse_channel}')]
        precursor_translated_pulse = df.loc[:, ('Precursor.Translated', f'{self.pulse_channel}')]
        precursor_quantity_pulse = df.loc[:, ('Precursor.Quantity', f'{self.pulse_channel}')]
        passes_filter_pulse = df.loc[:, ('filter_passed', f'{self.pulse_channel}')]
        
        # Combine into a new DataFrame
        combined_df = pd.DataFrame({
            'ms1_translated_pulse': ms1_translated_pulse,
            'precursor_translated_pulse': precursor_translated_pulse,
            'precursor_quantity_pulse' : precursor_quantity_pulse,
            'filter_passed_pulse' : passes_filter_pulse,
            
            'ms1_translated_L': ms1_translated_L,
            'precursor_translated_L': precursor_translated_L,
            'precursor_quantity_L' : precursor_quantity_L,
            'filter_passed_L' : passes_filter_L,
            
        })
        combined_df = combined_df.reset_index()
        return combined_df

    def format_data_dia_sis(self, df):
        # channel H
        ms1_translated_H = df.loc[:, ('Ms1.Translated', 'H')]
        precursor_translated_H = df.loc[:, ('Precursor.Translated', 'H')]
        precursor_quantity_H = df.loc[:, ('Precursor.Quantity', 'H')]
        passes_filter_H = df.loc[:, ('filter_passed', 'H')]
        
        # channel L
        ms1_translated_L = df.loc[:, ('Ms1.Translated', 'L')]
        precursor_translated_L = df.loc[:, ('Precursor.Translated', 'L')]
        precursor_quantity_L = df.loc[:, ('Precursor.Quantity', 'L')]
        passes_filter_L = df.loc[:, ('filter_passed', 'L')]
        
        # Combine into a new DataFrame
        combined_df = pd.DataFrame({
            'ms1_translated_H': ms1_translated_H,
            'precursor_translated_H': precursor_translated_H,
            'precursor_quantity_H' : precursor_quantity_H,
            'filter_passed_H' : passes_filter_H,
            
            'ms1_translated_L': ms1_translated_L,
            'precursor_translated_L': precursor_translated_L,
            'precursor_quantity_L' : precursor_quantity_L,
            'filter_passed_L' : passes_filter_L
        })
        combined_df = combined_df.reset_index()
        return combined_df

    def format_data_dynamic_dia_sis(self, df):
        # format new dataframes for future work (3 channels)
        # channel M
        ms1_translated_M = df.loc[:, ('Ms1.Translated', 'M')]
        precursor_translated_M = df.loc[:, ('Precursor.Translated', 'M')]
        precursor_quantity_M = df.loc[:, ('Precursor.Quantity', 'M')]
        passes_filter_M = df.loc[:, ('filter_passed', 'M')]
        
        # channel L
        ms1_translated_L = df.loc[:, ('Ms1.Translated', 'L')]
        precursor_translated_L = df.loc[:, ('Precursor.Translated', 'L')]
        precursor_quantity_L = df.loc[:, ('Precursor.Quantity', 'L')]
        passes_filter_L = df.loc[:, ('filter_passed', 'L')]
        
        # channel H
        ms1_translated_H = df.loc[:, ('Ms1.Translated', 'H')]
        precursor_translated_H = df.loc[:, ('Precursor.Translated', 'H')]
        precursor_quantity_H = df.loc[:, ('Precursor.Quantity', 'H')]
        passes_filter_H = df.loc[:, ('filter_passed', 'H')]
        
        # Combine into a new DataFrame
        combined_df = pd.DataFrame({
            'ms1_translated_H': ms1_translated_H,
            'precursor_translated_H': precursor_translated_H,
            'precursor_quantity_H' : precursor_quantity_H,
            'filter_passed_H' : passes_filter_H,
            
            'ms1_translated_M': ms1_translated_M,
            'precursor_translated_M': precursor_translated_M,
            'precursor_quantity_M' : precursor_quantity_M,
            'filter_passed_M' : passes_filter_M,
            
            'ms1_translated_L': ms1_translated_L,
            'precursor_translated_L': precursor_translated_L,
            'precursor_quantity_L' : precursor_quantity_L,
            'filter_passed_L' : passes_filter_L,
            
        })
        combined_df = combined_df.reset_index()
        return combined_df


