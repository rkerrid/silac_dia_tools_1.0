# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 18:27:52 2024

@author: robbi
"""

import pandas as pd
import operator


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
        file_path = f"{self.path}report.tsv"

        for count, chunk in enumerate(pd.read_table(file_path, sep="\t", chunksize=self.chunk_size), start=1):
            chunk['Genes'] = chunk['Genes'].fillna('')
            chunk['Protein.Group'] = chunk['Protein.Group'].str.cat(chunk['Genes'], sep='-')
            
            chunk = self.add_label_col(chunk)
            chunk = self.add_precursor_col(chunk)
            chunk = self.remove_cols(chunk)
            # print(chunk.columns.values.tolist())
        
            chunks.append(chunk)
            
            # if self.update:
            #     print(f'Chunk {count} processed')
            # if count == 10:
                # break
        df = pd.concat(chunks, ignore_index=True)
        print('Finished import')
        return df
    
    def add_label_col(self, chunk):
        chunk['Label'] = chunk['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]
        return chunk
    
    def add_precursor_col(self,chunk):
        chunk['Precursor'] = chunk['Stripped.Sequence'].astype(str) + chunk['Precursor.Charge'].astype(str)
        return chunk
    
    def remove_cols(self, chunk):
        cols = ['Run', 'Protein.Group', 'Precursor', 'Label',  'Precursor.Quantity','Ms1.Translated','Precursor.Translated'] + self.filter_cols
        chunk = chunk[cols]
        return chunk
    
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
