# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 18:20:58 2024

@author: robbi
"""

import pandas as pd
import numpy as np


class SilacFormatter:
    def __init__(self, path, filter_cols):
        self.path = path
        self.filter_cols = filter_cols

    def format_silac_channels(self, report):
        print('Beginning formatting SILAC channels')
        parsed_df = self._parse_data_for_channel_info(report)
        combined_precursors = self._combine_modified_precursors(parsed_df)
        stacked_intensities = self._stack_intensities(combined_precursors)
        print('Finished formatting SILAC channels')
        return stacked_intensities

    def _parse_data_for_channel_info(self, report):
        print('Parsing data for SILAC intensities')
        report['Label'] = report['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]
        report['Precursor'] = report['Stripped.Sequence'].astype(str) + report['Precursor.Charge'].astype(str)

        # Create Ms1.Translated and Precursor.Translated dataframes
        ms1_df, precursor_df = self._create_intensity_dfs(report)
        
        # Concatenate and return the parsed dataframe
        parsed_df = pd.concat([ms1_df, precursor_df], ignore_index=True)
        return parsed_df[['Run', 'Protein.Group', 'Precursor.Id', 'Precursor', 'Label', 'intensity', 'quantity type', 'Precursor.Quantity'] + self.filter_cols]

    def _create_intensity_dfs(self, report):
        ms1_df = report.copy()
        ms1_df['intensity'] = report['Ms1.Translated']
        ms1_df['quantity type'] = 'Ms1.Translated'
    
        precursor_df = report.copy()
        precursor_df['intensity'] = report['Precursor.Translated']
        precursor_df['quantity type'] = 'Precursor.Translated'
        return ms1_df, precursor_df

    def _combine_modified_precursors(self, parsed_df):
        print('Combining modified precursors')
        # Define aggregation functions for columns
        agg_functions = {key: 'first' for key in self.filter_cols}
        agg_functions['intensity'] = 'first'

        # Aggregate and pivot data
        pivoted_df = self._pivot_data(parsed_df, agg_functions)
        combined_df = self._merge_pivoted_data(parsed_df, pivoted_df)
        return combined_df.drop_duplicates()

    def _pivot_data(self, parsed_df, agg_functions):
        agg_df = parsed_df.groupby(['Run', 'Protein.Group', 'Precursor', 'quantity type', 'Label']).agg(agg_functions).reset_index()
        return agg_df.pivot_table(index=['Run', 'Protein.Group', 'Precursor', 'quantity type'], columns='Label', values='intensity', fill_value=0).reset_index()

    def _merge_pivoted_data(self, parsed_df, pivoted_df):
        # Rename columns and merge with original df
        pivoted_df.columns.name = None
        pivoted_df.rename(columns={'H': 'H intensity', 'M': 'M intensity', 'L': 'L intensity'}, inplace=True)
        merged_df = pd.merge(parsed_df.drop(columns=['Label', 'intensity']), pivoted_df, on=['Run', 'Protein.Group', 'Precursor', 'quantity type'])
        return merged_df.sort_values(self.filter_cols, ascending=False).drop_duplicates(subset=['Run', 'Protein.Group', 'Precursor', 'quantity type'])

    def _stack_intensities(self, combined_precursors):
        print("Stacking intensities")
        df = combined_precursors
        self._ensure_intensity_columns(df)
        self._calculate_ratios(df)
        df = self._drop_nan_rows(df)
        return self._select_final_columns(df)
    
    def _ensure_intensity_columns(self, df):
        for col in ['H intensity', 'M intensity', 'L intensity']:
            if col not in df.columns:
                df[col] = 0
    
    def _calculate_ratios(self, df):
        df['Precursor.Quantity'] = df['H intensity'] + df['M intensity'] + df['L intensity']
        for label in ['L', 'M', 'H']:
            ratio_column = f'{label} to stack ratio'
            df[ratio_column] = df[f'{label} intensity'] / df['Precursor.Quantity']
    
    def _drop_nan_rows(self, df):
        required_columns = [
            'Precursor.Quantity', 'H intensity', 'L intensity', 'M intensity', 
            'L to stack ratio', 'M to stack ratio', 'H to stack ratio'
        ]
        return df.dropna(subset=required_columns)
    
    def _select_final_columns(self, df):
        columns = ['Run', 'Protein.Group', 'Precursor.Id', 'Precursor.Quantity', 'quantity type', 
                   'H intensity', 'M intensity', 'L intensity', 'H to stack ratio', 'M to stack ratio', 'L to stack ratio'] + self.filter_cols
        return df[columns].rename(columns={'Precursor': 'Precursor.Id'})
    
