# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 18:40:32 2024

@author: robbi
"""

import pandas as pd
import numpy as np
from tqdm import tqdm

class PrecursorRollup:
    def __init__(self, path):
        self.path = path

    @staticmethod
    def _select_ms1_translated(group):
        return group[group['quantity type'] == 'Ms1.Translated']

    def calculate_protein_level_ratios(self, df):
        print("Calculating ratios from precursor information")
        protein_precursors = df.groupby(['Run', 'Protein.Group'])
        
        protein_data = []
        protein_count, protein_missed = 0, 0

        for name, group in tqdm(protein_precursors, desc="Processing proteins"):
            ms1_group = self._select_ms1_translated(group)
            if len(ms1_group) > 2:
                median_ratios = self._calculate_median_log_ratios(ms1_group)
                total_intensity = np.sum(ms1_group['Precursor.Quantity'])
                new_row = self._create_protein_row(group, median_ratios, total_intensity)
                protein_data.append(new_row)
                protein_count += 1
            else:
                protein_missed += 1

        self._print_protein_counts(protein_count, protein_missed, len(protein_precursors))
        return pd.DataFrame(protein_data)

    def _calculate_median_log_ratios(self, group):
        median_log2_ratios = np.median(np.log2(group[['L to stack ratio', 'M to stack ratio', 'H to stack ratio']]), axis=0)
        return np.exp2(median_log2_ratios)

    def _create_protein_row(self, group, median_ratios, total_intensity):
        return {
            'Run': group['Run'].iloc[0],
            'Protein.Group': group['Protein.Group'].iloc[0],
            'Total intensity': total_intensity,
            'L to stack ratio': median_ratios[0],
            'M to stack ratio': median_ratios[1],
            'H to stack ratio': median_ratios[2],
            'L intensity': median_ratios[0] * total_intensity,
            'M intensity': median_ratios[1] * total_intensity,
            'H intensity': median_ratios[2] * total_intensity
        }

    def _print_protein_counts(self, protein_count, protein_missed, total):
        print(f'Total proteins counted: {protein_count}')
        print(f'Total sets of precursors that didn\'t meet minimum unique precursor requirements: {protein_missed} out of {total}')

# Example usage:
# precursor_rollup = PrecursorRollup(path='your_path')
# df = pd.read_csv('path_to_dataframe.csv')
# protein_ratios = precursor_rollup.calculate_protein_level_ratios(df)
