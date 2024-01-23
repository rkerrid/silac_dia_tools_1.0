# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 18:27:52 2024

@author: robbi
"""

import pandas as pd
import numpy as np
import operator

from icecream import ic

class Preprocessor:
    def __init__(self, path, params, filter_cols, meta_data=None):
        self.path = path
        self.meta_data = meta_data
        self.contains_metadata = meta_data is not None
        self.params = params
        self.chunk_size = 10000
        self.update = True
        self.filter_cols = filter_cols 

    def import_report(self):
        print('Beginning import .tsv')
        chunks = []
        filtered_out = []
        file_path = f"{self.path}report.tsv"
        
        for count, chunk in enumerate(pd.read_table(file_path, sep="\t", chunksize=self.chunk_size), start=1):
            
            # reduce data size by subsetting report.tsv based on metadata, and removing columns not needed for further analysis
            # in the following loop we also annotate the silac chanels and append genes to Protein.Groups for downstream useage
            if self.meta_data is not None:
                chunk = self.subset_based_on_metadata(chunk)
            chunk['Genes'] = chunk['Genes'].fillna('')
            chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
            chunk = self.add_label_col(chunk)
            chunk = self.remove_cols(chunk)
            # annotate df with SILAC chanel then apply strict filters to H by droping the precursor, or adding NaN for L and M channels if they dont pass loose filters
            chunk, chunk_filtered_out = self.filter_spikein_strict(chunk, "H")
            # chunk = self.apply_nan_by_loose_filtering(chunk,"L")
            # chunk = self.apply_nan_by_loose_filtering(chunk,"M")
            
            chunks.append(chunk)
            filtered_out.append(chunk_filtered_out)
            
            if self.update:
                print(f'Chunk {count} processed')
            # if count == 1:
            #     break
            
        df = pd.concat(chunks, ignore_index=True)
        filtered_out_df = pd.concat(filtered_out, ignore_index=True)
        print('Finished import')
        return df, filtered_out_df
    
    def subset_based_on_metadata(self, chunk):
        # need to complete this method
        return chunk
    
    def add_label_col(self, chunk):
        chunk['Label'] = chunk['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]
        return chunk
    
    
    def remove_cols(self, chunk):
        cols = ['Run', 'Protein.Group', 'Precursor.Id', 'Label',  'Precursor.Quantity','Ms1.Translated','Precursor.Translated'] + self.filter_cols 
        chunk = chunk[cols]
        return chunk
    
    def filter_spikein_strict(self, chunk, label):
        # Initialize boolean mask the size of the df to be switched to true or false when looping through columns and conditions in params
        filtering_condition = pd.Series([True] * len(chunk), index=chunk.index)
        chunk_filtered_out = chunk.copy()
        ops = {
            "==": operator.eq, "<": operator.lt, "<=": operator.le,
            ">": operator.gt, ">=": operator.ge
        }

        if label in chunk['Label'].values:
            for column, condition in self.params['apply_strict_filters'].items():
                op = ops[condition['op']]
                # &= means if the row hasn't passed a previous filter it will remain False
                filtering_condition &= op(chunk[column], condition['value'])
            # ic(filtering_condition)
            # chunk_copy = chunk.copy(deep=True)
            filtered_chunk = chunk[filtering_condition]
            chunk_filtered_out = chunk[~filtering_condition]
        else:
            # If the label is not present, return the whole chunk and an empty DataFrame
            filtered_chunk = chunk
            chunk_filtered_out = pd.DataFrame(columns=chunk.columns)

        return filtered_chunk, chunk_filtered_out
    
    

    def apply_nan_by_loose_filtering(self, chunk, label):
        filtering_condition = pd.Series([True] * len(chunk), index=chunk.index)
        ops = {
            "==": operator.eq, "<": operator.lt, "<=": operator.le,
            ">": operator.gt, ">=": operator.ge
        }
        filtered_chunk = chunk.copy()
    
        if label in chunk['Label'].values:
            for column, condition in self.params['apply_loose_filters'].items():
                op = ops[condition['op']]
                # Update the condition for each column
                filtering_condition &= op(chunk[column], condition['value'])
    
            nan_cols = ['Precursor.Translated', 'Precursor.Quantity', 'Ms1.Translated']
            for col in nan_cols:
                if col in chunk.columns:
                    # Set NaN where the condition is False
                    filtered_chunk.loc[~filtering_condition, col] = np.nan
    
        return filtered_chunk


    # def apply_nan_by_loose_filtering(self, chunk, label):
    #     filtering_condition = pd.Series([False] * len(chunk), index=chunk.index)
    #     ops = {
    #         "==": operator.eq, "<": operator.lt, "<=": operator.le,
    #         ">": operator.gt, ">=": operator.ge
    #     }
    #     filtered_chunk = chunk.copy()

    #     if label in chunk['Label']:
    #         for column, condition in self.params['apply_strict_filters'].items():
                
    #             op = ops[condition['op']]
    #             filtering_condition &= op(chunk[column], condition['value'])
                
    #             nan_cols = ['Precursor.Translated', 'Precursor.Quantity', 'Ms1.Translated']
    #             for col in nan_cols:
    #                 if col in chunk.columns:
    #                     filtered_chunk.loc[filtering_condition, col] = np.nan
    
    #     return filtered_chunk
    
    
    
    
    
    
    
    def filter_formatted(self, formatted_precursors):
        print('Begin filtering formatted precursors')
        if self.contains_metadata:
            formatted_precursors = self._relable_run(formatted_precursors)
        cleaned_precursors, contaminants = self._remove_contaminants(formatted_precursors)
        filtered_precursors, filtered_out = self._apply_filters(cleaned_precursors)
        return filtered_precursors, contaminants, filtered_out

    def _remove_contaminants(self, chunk):
        contams_mask = chunk['Protein.Group'].str.contains('Cont_', case=False, na=False)
        self._validate_boolean_mask(contams_mask)
        
        contaminants = chunk[contams_mask]
        cleaned_chunk = chunk[~contams_mask]
        return cleaned_chunk, contaminants

    def _apply_filters(self, chunk):
        filtering_condition = pd.Series([True] * len(chunk), index=chunk.index)
        ops = {
            "==": operator.eq, "<": operator.lt, "<=": operator.le,
            ">": operator.gt, ">=": operator.ge
        }

        for column, condition in self.params['apply_filters'].items():
            op = ops[condition['op']]
            filtering_condition &= op(chunk[column], condition['value'])

        filtered_chunk = chunk[filtering_condition]
        chunk_filtered_out = chunk[~filtering_condition]
        return filtered_chunk, chunk_filtered_out

    def _relable_run(self, chunk):
        run_to_sample = dict(zip(self.meta_data['Run'], self.meta_data['Sample']))
        chunk['Run'] = chunk['Run'].map(run_to_sample)
        if chunk['Run'].isna().any():
            raise ValueError("Some Run values in report.tsv not found in metadata.")
        return chunk

    def _validate_boolean_mask(self, mask):
        if not all(isinstance(x, bool) for x in mask):
            invalid_values = mask[~mask.isin([True, False])]
            print(f"Non-boolean values in mask: {invalid_values}")

    def drop_cols(self, chunk, filter_cols=[]):
        chunk['Genes'] = chunk['Genes'].fillna('')
        chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
        cols_to_keep = ['Run', 'Protein.Group', 'Stripped.Sequence', 'Precursor.Id', 'Precursor.Charge',
                        'Precursor.Quantity', 'Precursor.Translated', 'Ms1.Translated'] + filter_cols
        return chunk[cols_to_keep]
