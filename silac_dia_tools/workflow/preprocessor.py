# # -*- coding: utf-8 -*-
# """
# Created on Tue Sep  3 10:14:26 2024

# @author: robbi
# """

# import pandas as pd
# import numpy as np
# import operator
# import time 

# from tqdm import tqdm
# import os
# from icecream import ic


# class Preprocessor:
#     def __init__(self, path, method, pulse_channel=None, meta_data=None):
#         self.path = path
#         self.meta_data = meta_data
#         self.chunk_size = 180000
     
#         self.pulse_channel = pulse_channel
#         self.method = method
        
#     def import_report(self):
#         print('Beginning import report.tsv')
#         start_time = time.time()
                
#         chunks = []
#         contaminants = []
#         file_path = f"{self.path}report.tsv"
#         count = 1
#         # Estimate rows in file size
#         file_size_bytes = os.path.getsize(file_path)
#         average_row_size_bytes = 1000  # This is an example; you'll need to adjust this based on your data
#         # Estimate the number of rows
#         estimated_rows = file_size_bytes / average_row_size_bytes
#         total_chunks = estimated_rows/self.chunk_size
#         for chunk in tqdm(pd.read_table(file_path, sep="\t", chunksize=self.chunk_size), 
#                       total=total_chunks, desc='Estimated loading of report.tsv based on file size'):
         
#             if self.meta_data is not None:
#                 chunk = self.subset_based_on_metadata(chunk)
#                 chunk = self.relabel_run(chunk)
                
#             chunk = self.add_label_col(chunk)
#             chunk = self.add_passes_filter_col(chunk)
#             chunk = self.drop_cols(chunk)
#             chunk = self.pivot_data(chunk)
            
#             chunk = self.rename_cols(chunk, self.method, self.pulse_channel)
           
#             chunk, contam_chunk = self.remove_contaminants(chunk)
            
#             chunks.append(chunk)
#             contaminants.append(contam_chunk)
            
#         # append chunks to respective dfs and return  
#         filtered_report = pd.concat(chunks, ignore_index=True)
#         contaminants_df = pd.concat(contaminants, ignore_index=True)
#         print('Finished import')
#         end_time = time.time()
#         print(f"Time taken for import: {end_time - start_time} seconds")
        
#         return filtered_report, contaminants_df


#     def subset_based_on_metadata(self, df):       
#         filtered_df = df[df['Run'].isin(self.meta_data['Run'])]
#         return filtered_df
    
#     def relabel_run(self, df):
#         run_to_sample = dict(zip(self.meta_data['Run'], self.meta_data['Sample']))
    
#         # Apply the mapping to df2['Run'] and raise an error if a 'Run' value doesn't exist in df1
#         df['Run'] = df['Run'].map(run_to_sample)
#         if df['Run'].isna().any():
#             raise ValueError("Some Run values in the report.tsv are not found in the metadata, please ensure metadata is correct.")
            
#         return df

#     def add_label_col(self, df):
#         # Extract the label and add it as a new column
#         df['Label'] = df['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]
        
#         # Remove the '(SILAC-K|R-([HML]))' part from the 'Precursor.Id' string
#         df['Precursor.Id'] = df['Precursor.Id'].str.replace(r'\(SILAC-(K|R)-[HML]\)', '', regex=True)
    
#         return df

#     def add_passes_filter_col(self, df):
        
#         # if data does not pass the following filters set to False
#         df['filter_passed'] = (df["Global.PG.Q.Value"] < 0.01) & (df["Precursor.Charge"] > 1) & (df["Channel.Q.Value"] < 0.03)
    
#         return df
    
#     def drop_cols(self, df):
#         # what cols to keep for future workflow
#         cols = ['Run',
#                  'Protein.Group',
#                  'Protein.Ids',
#                  'Protein.Names',
#                  'Genes',
#                  'Precursor.Id',
#                  'Precursor.Quantity',
#                  'Precursor.Normalised',
#                  'Precursor.Translated',
#                  'Ms1.Translated',
#                  'Label',
#                  'filter_passed']
        
#         # drop all other cols
#         df = df[cols]
#         return df


#     def remove_contaminants(self, df):
#         #chunk_copy = chunk.copy(deep=True)
#         contams_mask = df['protein_group'].str.contains('Cont_', case=False, na=False)
#         df_filtered = df.loc[~contams_mask].reset_index(drop=True)
#         contams = df.loc[contams_mask].reset_index(drop=True)
   
#         return df_filtered, contams
    
    
#     #pivot data
#     def pivot_data(self, df):
#         index_cols = ['Run', 'Protein.Group','Protein.Ids','Protein.Names', 'Genes', 'Precursor.Id'] #, 'Protein.Ids', 'Protein.Names', 'Genes', 'Precursor.Id', 'Passes_filter'
#         df_p = df.pivot_table(index=index_cols, 
#                               columns='Label', 
#                               values = ['Ms1.Translated', 'Precursor.Quantity', 'Precursor.Translated','filter_passed'],
#                               aggfunc=lambda x: x)
#         df_p['filter_passed'] = df_p['filter_passed'].applymap(lambda x: x == 1.0)
        
#         # Flatten the multi-level columns
#         df_p.columns = [f'{col[0]}_{col[1]}' for col in df_p.columns]
        
#         # Reset index to get a clean DataFrame
#         df_p = df_p.reset_index()
#         return df_p
    
#     def rename_cols(self, df, channel, method):
#         df.rename(columns={'Protein.Group':'protein_group','Protein.Ids':'protein_ids', 'Protein.Names':'protein_names', 'Genes':'genes', 'Precursor.Id': 'precursor_id'}, inplace=True)
#         if method == 'dynamic_silac_dia': 
#             df.rename(columns={'Ms1.Translated_L':'ms1_translated_L',
#                                f'Ms1.Translated_{channel}':'ms1_translated_pulse', 
#                                'Precursor.Translated_L':'precursor_translated_L',
#                                f'Precursor.Translated_{channel}':'precursor_translated_pulse',
#                                'Precursor.Quantity_L':'precursor_quantity_L',
#                                f'Precursor.Quantity_{channel}':'precursor_quantity_pulse',
#                                f'filter_passed_{channel}':'filter_passed_pulse'}, inplace=True)
        
#         elif method == 'dynamic_dia_sis':
#             df.rename(columns={'Ms1.Translated_L':'ms1_translated_L',
#                                f'Ms1.Translated_{channel}':'ms1_translated_pulse', 
#                                'Precursor.Translated_L':'precursor_translated_L',
#                                f'Precursor.Translated_{channel}':'precursor_translated_pulse',
#                                'Precursor.Quantity_L':'precursor_quantity_L',
#                                f'Precursor.Quantity_{channel}':'precursor_quantity_pulse',
#                                f'filter_passed_{channel}':'filter_passed_pulse',
#                                'Ms1.Translated_H':'ms1_translated_H',
#                                'Precursor.Translated_H':'precursor_translated_H',
#                                'Precursor.Quantity_H':'precursor_quantity_H'}, inplace=True)
            
#         elif method == 'dia_sis':
#             df.rename(columns={'Ms1.Translated_L':'ms1_translated_L',
#                                'Ms1.Translated_H':'ms1_translated_H', 
#                                'Precursor.Translated_L':'precursor_translated_L',
#                                'Precursor.Translated_H':'precursor_translated_H',
#                                'Precursor.Quantity_L':'precursor_quantity_L',
#                                'Precursor.Quantity_H':'precursor_quantity_H',
#                                'filter_passed_H':'filter_passed_H'}, inplace=True)
#         return df
    
     
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
from concurrent.futures import ProcessPoolExecutor


# class Preprocessor:
#     def __init__(self, path, method, pulse_channel=None, meta_data=None):
#         self.path = path
#         self.meta_data = meta_data
#         self.chunk_size = 180000
     
#         self.pulse_channel = pulse_channel
#         self.method = method
        
#     def import_report(self):
#         print('Beginning import report.tsv')
#         start_time = time.time()
                
#         chunks = []
#         contaminants = []
#         file_path = f"{self.path}report.tsv"
#         count = 1
#         # Estimate rows in file size
#         file_size_bytes = os.path.getsize(file_path)
#         average_row_size_bytes = 1000  # This is an example; you'll need to adjust this based on your data
#         # Estimate the number of rows
#         estimated_rows = file_size_bytes / average_row_size_bytes
#         total_chunks = estimated_rows/self.chunk_size
#         for chunk in tqdm(pd.read_table(file_path, sep="\t", chunksize=self.chunk_size), 
#                       total=total_chunks, desc='Estimated loading of report.tsv based on file size'):
         
#             if self.meta_data is not None:
#                 chunk = self.subset_based_on_metadata(chunk)
#                 chunk = self.relabel_run(chunk)
                
#             chunk = self.add_label_col(chunk)
#             chunk = self.add_passes_filter_col(chunk)
#             chunk = self.drop_cols(chunk)
#             chunk = self.pivot_data(chunk)
            
#             chunk = self.rename_cols(chunk, self.method, self.pulse_channel)
           
#             chunk, contam_chunk = self.remove_contaminants(chunk)
            
#             chunks.append(chunk)
#             contaminants.append(contam_chunk)
            
#         # append chunks to respective dfs and return  
#         filtered_report = pd.concat(chunks, ignore_index=True)
#         contaminants_df = pd.concat(contaminants, ignore_index=True)
#         print('Finished import')
#         end_time = time.time()
#         print(f"Time taken for import: {end_time - start_time} seconds")
        
#         return filtered_report, contaminants_df

class Preprocessor:
    def __init__(self, path, method, pulse_channel=None, meta_data=None):
        self.path = path
        self.meta_data = meta_data
        self.chunk_size = 1000000  # Adjusted chunk size for better performance
        self.pulse_channel = pulse_channel
        self.method = method
    
    def preprocess(self):
        filtered_report, contaminants_df = self.import_report()
        filtered_report.to_csv(f'{self.path}to_format.csv', sep=',')
        print('Reformating')
        start_time = time.time()
        
        filtered_report = self.reformat_table(filtered_report, self.method, self.pulse_channel)
        ic(filtered_report)
        # filtered_report = self.pivot_data(filtered_report)
        # filtered_report = self.rename_cols(filtered_report, self.method, self.pulse_channel)
        print('Finished reformating')
        end_time = time.time()
        print(f"Time taken for reformating: {end_time - start_time} seconds")
        return filtered_report, contaminants_df
    
    def reformat_table(self, df, method, pulse_channel):
        df = df.rename(columns={'Protein.Group':'protein_group','Protein.Ids':'protein_ids', 'Protein.Names':'protein_names', 'Genes':'genes', 'Precursor.Id': 'precursor_id'})
        index_cols = ['Run', 'protein_group', 'protein_ids', 'protein_names', 'genes', 'precursor_id']
        if method == 'dynamic_silac_dia':         
            df_light = df[df['Label']=='L']
            df_pulse = df[df['Label']==pulse_channel]
            df_light = df_light.drop(['Label'], axis=1)
            df_pulse = df_pulse.drop(['Label'], axis=1)
     
            df_light = df_light.rename(columns={'Precursor.Quantity':'precursor_quantity_L','Precursor.Translated':'precursor_translated_L','Ms1.Translated':'ms1_translated_L','filter_passed':'filter_passed_L'})
            df_pulse = df_pulse.rename(columns={'Precursor.Quantity':'precursor_quantity_pulse','Precursor.Translated':'precursor_translated_pulse','Ms1.Translated':'ms1_translated_pulse','filter_passed':'filter_passed_pulse'})
           
            df = pd.merge(df_light, df_pulse,on=index_cols, how='outer')
           
            return df
        
        elif method == 'dynamic_dia_sis':
            df_L = df[df['Label']=='L']
            df_M = df[df['Label']=='M']
            df_H = df[df['Label']=='H']
            
            df_L = df_L.drop(['Label'], axis=1)
            df_M = df_M.drop(['Label'], axis=1)
            df_H = df_H.drop(['Label'], axis=1)

            df_L = df_L.rename(columns={'Precursor.Quantity':'precursor_quantity_L','Precursor.Translated':'precursor_translated_L','Ms1.Translated':'ms1_translated_L','filter_passed':'filter_passed_L'})
            df_M = df_M.rename(columns={'Precursor.Quantity':'precursor_quantity_M','Precursor.Translated':'precursor_translated_M','Ms1.Translated':'ms1_translated_M','filter_passed':'filter_passed_M'})
            df_H = df_H.rename(columns={'Precursor.Quantity':'precursor_quantity_H','Precursor.Translated':'precursor_translated_H','Ms1.Translated':'ms1_translated_H','filter_passed':'filter_passed_H'})
            
            df = df_L.merge(df_M, on=index_cols, how='outer').merge(df_H, on=index_cols, how='outer')
            df['filter_passed_H'] = df['filter_passed_H'].fillna(False)
      
            return df
    
    def import_report(self):
        print('Beginning import of report.tsv')
        start_time = time.time()
        
        file_path = f"{self.path}report.tsv"
      
        file_size_bytes = os.path.getsize(file_path)
        average_row_size_bytes = 1100  # This is an example; you'll need to adjust this based on your data
        # Estimate the number of rows
        estimated_rows = file_size_bytes / average_row_size_bytes
        total_chunks = estimated_rows/self.chunk_size
        # Use ProcessPoolExecutor for parallel processing
        with ProcessPoolExecutor() as executor:
            futures = []
            with tqdm(total=total_chunks, desc="Processing file in chunks") as pbar:
                for chunk in pd.read_table(file_path, sep="\t", chunksize=self.chunk_size):
                    futures.append(executor.submit(self.process_chunk, chunk))
                    pbar.update(1) 
            
            # Gather results from futures
            results = [f.result() for f in futures]
            
            # Concatenate chunks into final DataFrame
            filtered_report = pd.concat([res[0] for res in results], ignore_index=True)
            contaminants_df = pd.concat([res[1] for res in results], ignore_index=True)

        print('Finished import')
        end_time = time.time()
        print(f"Time taken for import: {end_time - start_time} seconds")
        
        return filtered_report, contaminants_df

    def process_chunk(self, chunk):
        # Process the chunk in parallel
        if self.meta_data is not None:
            chunk = self.subset_based_on_metadata(chunk)
            chunk = self.relabel_run(chunk)
        
        chunk = self.add_label_col(chunk)
        chunk = self.add_passes_filter_col(chunk)
        chunk = self.drop_cols(chunk)
        
        chunk, contam_chunk = self.remove_contaminants(chunk)
        return chunk, contam_chunk
    
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
                  'Precursor.Translated',
                  'Ms1.Translated',
                  'Label',
                  'filter_passed']
        
        # drop all other cols
        df = df[cols]
        return df


    def remove_contaminants(self, df):
        #chunk_copy = chunk.copy(deep=True)
        contams_mask = df['Protein.Group'].str.contains('Cont_', case=False, na=False)
        df_filtered = df.loc[~contams_mask].reset_index(drop=True)
        contams = df.loc[contams_mask].reset_index(drop=True)
   
        return df_filtered, contams
    
    
    #pivot data
    def pivot_data(self, df):
        index_cols = ['Run', 'Protein.Group','Protein.Ids','Protein.Names', 'Genes', 'Precursor.Id'] #, 'Protein.Ids', 'Protein.Names', 'Genes', 'Precursor.Id', 'Passes_filter'
        df_p = df.pivot_table(index=index_cols, 
                              columns='Label', 
                              values = ['Ms1.Translated', 'Precursor.Quantity', 'Precursor.Translated','filter_passed'],
                              aggfunc=lambda x: x)
        df_p['filter_passed'] = df_p['filter_passed'].applymap(lambda x: x == 1.0)
        
        # Flatten the multi-level columns
        df_p.columns = [f'{col[0]}_{col[1]}' for col in df_p.columns]
        
        # Reset index to get a clean DataFrame
        df_p = df_p.reset_index()
        return df_p
    
    def rename_cols(self, df, channel, method):
        df.rename(columns={'Protein.Group':'protein_group','Protein.Ids':'protein_ids', 'Protein.Names':'protein_names', 'Genes':'genes', 'Precursor.Id': 'precursor_id'}, inplace=True)
        if method == 'dynamic_silac_dia': 
            df.rename(columns={'Ms1.Translated_L':'ms1_translated_L',
                                f'Ms1.Translated_{channel}':'ms1_translated_pulse', 
                                'Precursor.Translated_L':'precursor_translated_L',
                                f'Precursor.Translated_{channel}':'precursor_translated_pulse',
                                'Precursor.Quantity_L':'precursor_quantity_L',
                                f'Precursor.Quantity_{channel}':'precursor_quantity_pulse',
                                f'filter_passed_{channel}':'filter_passed_pulse'}, inplace=True)
        
        elif method == 'dynamic_dia_sis':
            df.rename(columns={'Ms1.Translated_L':'ms1_translated_L',
                                f'Ms1.Translated_{channel}':'ms1_translated_M', 
                                'Precursor.Translated_L':'precursor_translated_L',
                                f'Precursor.Translated_{channel}':'precursor_translated_M',
                                'Precursor.Quantity_L':'precursor_quantity_L',
                                f'Precursor.Quantity_{channel}':'precursor_quantity_M',
                                f'filter_passed_{channel}':'filter_passed_M',
                                'Ms1.Translated_H':'ms1_translated_H',
                                'Precursor.Translated_H':'precursor_translated_H',
                                'Precursor.Quantity_H':'precursor_quantity_H'}, inplace=True)
            
        elif method == 'dia_sis':
            df.rename(columns={'Ms1.Translated_L':'ms1_translated_L',
                                'Ms1.Translated_H':'ms1_translated_H', 
                                'Precursor.Translated_L':'precursor_translated_L',
                                'Precursor.Translated_H':'precursor_translated_H',
                                'Precursor.Quantity_L':'precursor_quantity_L',
                                'Precursor.Quantity_H':'precursor_quantity_H'}, inplace=True)
        return df
    
     
