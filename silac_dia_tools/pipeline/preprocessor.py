# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 11:43:50 2023

@author: rkerrid

Step 1: Module for filtering report.tsv output (DIA-NN version 1.8.1) with
SILAC settings as described in the README.md

Note: This script filters for contaminants by looking for the 'cont_' substring
in Protein.Groups so make sure your report.tsv is annotated in the same way or 
edit the remove_contaminants() funciton.

"""

import pandas as pd
import json
import os
import operator

from icecream import ic
ic.disable()


class Preprocessor:
    def __init__(self, path, params, meta_data=None):
        self.path = path
        self.meta_data = meta_data
        self.contains_metadata = False
        if self.meta_data is not None:
            self.contains_metadata = True
        self.params = params
        self.filter_cols = self.params['apply_filters'].keys()
        self.chunk_size = 10000
        self.update = True
        
    def import_report(self):
        print('Beggining import no filter')
        
        count = 1
        with open(f"{self.path}report.tsv", 'r', encoding='utf-8') as file:
            chunks = []
    
            for chunk in pd.read_table(file,sep="\t", chunksize=self.chunk_size):
                chunk = self.subset_import(chunk)
                chunk['Genes'] = chunk['Genes'].fillna('')
                chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
                chunks.append(chunk)
                
                # Update progress (optional)
                if self.update:
                    print('chunk ', count,' processed')
                count+=1
            
            # Concatenate all chunks into a DataFrames
  
            df = pd.concat(chunks, ignore_index=True)
            print('Finished import')
        print("save subset")
        df.to_csv("{self.path}report_subset.tsv", sep=',')
        return df
    
    #remove samples not in the metadata from the report.tsv
    def subset_import(self, chunk):
        
        chunk = chunk[chunk['Run'].isin(self.meta_data['Run'])]
        
        return chunk
    
    def drop_non_meta_samples(self, chunk, meta):
        filtered_chunk = chunk[chunk['Run'].isin(meta['Run'])]
        return filtered_chunk
        
    
    def filter_formatted(self, formatted_precursors):
        print('Begin filtering formatted precursors')
        if self.contains_metadata:
            precursors = self.relable_run(formatted_precursors)
        precursors, contam = self.remove_contaminants(precursors)
        precursors, filtered_out = self.apply_filters(precursors)
     
        return precursors, contam, filtered_out
    
    #Filtering
    def remove_contaminants(self, chunk): # is self needed?
        # Create a contaminants mask based on the cont_ string and make sure all values are boolean
        contams_mask = chunk['Protein.Group'].str.contains('Cont_', case=False, na=False)
        if not all(isinstance(x, bool) for x in contams_mask):
            print("contams_mask contains non-boolean values:", contams_mask[~contams_mask.isin([True, False])])
        
        contams_df = chunk[contams_mask]  # Dataframe with only contaminants
        cleaned_chunk = chunk[~contams_mask]  # Dataframe without contaminants
        return cleaned_chunk, contams_df
        
    def apply_filters(self, chunk):
        # Initialize operator dict
        ops = {
            "==": operator.eq,
            "<": operator.lt,
            "<=": operator.le,
            ">": operator.gt,
            ">=": operator.ge
        }

         # Create a boolean Series with all True values and explicitly set its index
        filtering_condition = pd.Series([True] * len(chunk), index=chunk.index)
        
        # Iterating over each filter condition in params['apply_filters']
        for column, condition in self.params['apply_filters'].items():
            op = ops[condition['op']]
            value = condition['value']
            
            # Updating filtering_condition by applying each condition
            filtering_condition &= op(chunk[column], value)

        # Filter chunk and return both filtered and filtered out dfs
        chunk_filtered = chunk[filtering_condition]
        chunk_filtered_out = chunk[~filtering_condition]

        return chunk_filtered, chunk_filtered_out

    def  drop_cols(self, chunk, filter_cols = []): # is self needed?
        chunk['Genes'] = chunk['Genes'].fillna('')
        chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
        cols_to_keep = [ 'Run',
                          'Protein.Group',
                          'Stripped.Sequence',
                          'Precursor.Id', 
                          'Precursor.Charge',
            
                          'Precursor.Quantity',
                          'Precursor.Translated',
                          'Ms1.Translated'
                          ] + filter_cols
        chunk = chunk[cols_to_keep]
        return chunk
    
    def relable_run(self, chunk):
        run_to_sample = dict(zip(self.meta_data['Run'], self.meta_data['Sample']))

        # Apply the mapping to df2['Run'] and raise an error if a 'Run' value doesn't exist in df1
        chunk['Run'] = chunk['Run'].map(run_to_sample)
        if chunk['Run'].isna().any():
            raise ValueError("Some Run values in the report.tsv are not found in the metadata, please ensure metadata is correct.")
            
        return chunk
 