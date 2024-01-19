# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 18:32:12 2024

@author: robbi
"""

import pandas as pd
from silac_dia_tools.pipeline.utils import dlfq_functions as dlfq
from silac_dia_tools.pipeline.utils import manage_directories


class IntensityCalculator:
    def __init__(self, path, contains_reference, pulse_channel):
        self.path = path
        self.contains_reference = contains_reference
        self.pulse_channel = pulse_channel

    def output_href(self, ratios):
        manage_directories.create_directory(self.path, 'protein intensities')
        print('Calculating href intensities')
        href_normalized = self._calculate_href_normalized(ratios)
        return self._save_intensity_csvs(href_normalized, 'href')

    def _calculate_href_normalized(self, ratios):
        h_ref = self._get_h_reference(ratios)
        merged_df = ratios.merge(h_ref, on='Protein.Group', how='inner')
        self._assign_href_intensities(merged_df)
        return merged_df

    def _get_h_reference(self, ratios):
        h_ref = ratios.groupby('Protein.Group')['H intensity'].median()
        return h_ref.reset_index().rename(columns={'H intensity': 'h_ref'})

    def _assign_href_intensities(self, merged_df):
        merged_df['H normalized total intensity'] = merged_df['h_ref'] / merged_df['H to stack ratio']
        for channel in ['H', 'M', 'L']:
            merged_df[f'{channel} intensity'] = merged_df['H normalized total intensity'] * merged_df[f'{channel} to stack ratio']
        self._calculate_total_and_nsp_intensities(merged_df)

    def output_unnorm(self, ratios):
        manage_directories.create_directory(self.path, 'protein intensities')
        print('Calculating unnormalized intensities')
        unnorm = self._calculate_unnormalized(ratios)
        return self._save_intensity_csvs(unnorm, 'unnorm')

    def _calculate_unnormalized(self, ratios):
        nsp_channel = f'{self.pulse_channel} intensity'
        ratios['Total intensity'] = ratios['L intensity'] + ratios[nsp_channel]
        ratios['NSP intensity'] = ratios[nsp_channel]
        return ratios

    def output_dlfq(self, ratios, silac_precursors):
        manage_directories.create_directory(self.path, 'protein intensities')
        print('Calculating dlfq intensities')
        dlfq_normalized = self._calculate_dlfq_normalized(ratios, silac_precursors)
        return self._save_intensity_csvs(dlfq_normalized, 'dlfq')

    def _calculate_dlfq_normalized(self, ratios, silac_precursors):
        self._prepare_and_run_directlfq(silac_precursors)
        lfq_df = self._load_directlfq_results()
        merged_df = ratios.merge(lfq_df, on=['Protein.Group', 'Run'], how='inner')
        self._assign_dlfq_intensities(merged_df)
        return merged_df

    def _prepare_and_run_directlfq(self, silac_precursors):
        silac_precursors_file = f'{self.path}preprocessing/silac_precursors_dlfq_in.tsv'
        dlfq_output_file = f'{self.path}preprocessing/dlfq_protein_intensities.tsv'
        silac_precursors.to_csv(silac_precursors_file, sep='\t')
        dlfq.run_lfq(silac_precursors_file, file=dlfq_output_file, num_cores=1)

    def _load_directlfq_results(self):
        lfq_df = pd.read_csv(f'{self.path}preprocessing/dlfq_protein_intensities.tsv', sep='\t')
        return lfq_df.melt(id_vars=['protein'], var_name='Run', value_name='Intensity').rename(columns={'protein': 'Protein.Group'})

    def _assign_dlfq_intensities(self, merged_df):
        nsp_ratio = f'{self.pulse_channel} to stack ratio'
        for channel in ['L', 'Total', 'NSP']:
            if channel == 'Total':
                merged_df['Total intensity'] = (merged_df['L to stack ratio'] + merged_df[nsp_ratio]) * merged_df['Intensity']
            else:
                merged_df[f'{channel} intensity'] = merged_df[f'{channel[0]} to stack ratio'] * merged_df['Intensity']
    
    def _calculate_total_and_nsp_intensities(self, df):
        df['Total intensity'] = df['L intensity'] + df['M intensity']
        df['NSP intensity'] = df['M intensity']
    
    def _save_intensity_csvs(self, df, method):
        subsets = {
            'light': 'L intensity',
            'nsp': 'NSP intensity',
            'total': 'Total intensity'
        }
        if self.contains_reference or method != 'unnorm':
            subsets['reference'] = 'H intensity'
    
        results = {}
        for key, column in subsets.items():
            pivoted_df = df.pivot(index='Protein.Group', columns='Run', values=column)
            file_path = f'{self.path}protein intensities/{key}_{method}.csv'
            pivoted_df.to_csv(file_path, sep=',')
            results[key] = pivoted_df
    
        print(f'Saved {method} normalized protein intensities')
        return results['nsp'], results['total'], results.get('light')
