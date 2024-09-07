# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 08:55:46 2024

@author: robbi
"""


import unittest
from silac_dia_tools.workflow.dynamic_silac_dia import DynamicSilac
import pandas as pd
import pandas.testing as pdt
from icecream import ic

class TestDynamicSilac(unittest.TestCase):
 
    def test_initialize_class(self):
        report = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_reformatted_no_spike.csv', sep=',')
        precursor_rollup = DynamicSilac('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/', report)
    
    def test_filter_data(self):
        report = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_reformatted_no_spike.csv', sep=',')
        precursor_rollup = DynamicSilac('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/', report)
        
        report[['filter_passed_L','filter_passed_pulse']] = report[['filter_passed_L','filter_passed_pulse']].astype('bool')
        df = precursor_rollup.filter_data(report)
       
        # df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/filtered_no_spike.csv', index=False)
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/filtered_no_spike.csv', sep=',')
      
        
        self.assertEqual(len(df), len(df_compare), "DataFrames should have the same number of rows")
        
    def test_calculate_ratios(self):
        report = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_reformatted_no_spike.csv', sep=',')
        precursor_rollup = DynamicSilac('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/', report)
        report[['filter_passed_L','filter_passed_pulse']] = report[['filter_passed_L','filter_passed_pulse']].astype('bool')
        
        df = precursor_rollup.filter_data(report)
        df = precursor_rollup.calculate_precursor_ratios(df)
        df = df.reset_index(drop=True)
       
        
        # df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/precursor_ratios_no_spike.csv', index=False)
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/precursor_ratios_no_spike.csv', sep=',')
        
        df_compare[['filter_passed_L','filter_passed_pulse']] = df_compare[['filter_passed_L','filter_passed_pulse']].astype('bool')
        
        pdt.assert_frame_equal(df, df_compare)
    
    def test_compute_protein_level_ratios(self):
        df = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/precursor_ratios_no_spike.csv', sep=',')
        precursor_rollup = DynamicSilac('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/', df)
        
        df = precursor_rollup.compute_protein_level_ratios(df)
        # df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_group_ratios.csv', sep=',', index=False)
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_group_ratios.csv', sep=',')
        
        pdt.assert_frame_equal(df, df_compare)
        
    def test_perform_lfq(self):
        df = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/precursor_ratios_no_spike.csv', sep=',')
        df[['filter_passed_L','filter_passed_pulse']] = df[['filter_passed_L','filter_passed_pulse']].astype('bool')
        
        precursor_rollup = DynamicSilac('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/', df)
        
        lfq_df = precursor_rollup.perform_lfq(df)
        lfq_df = lfq_df.reset_index(drop=True)
        
        # lfq_df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_group_intensities.csv', sep=',', index=False)
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_group_intensities.csv', sep=',')
        
        pdt.assert_frame_equal(lfq_df, df_compare)
        
    def test_merge_data(self):
        ratios = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_group_ratios.csv', sep=',')
        intensities = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_group_intensities.csv', sep=',')
        precursor_rollup = DynamicSilac('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/', ratios)
        
        df = precursor_rollup.merge_data(ratios, intensities)
        df = df.reset_index(drop=True)
        
        # df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_group_merged.csv', sep=',', index=False)
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_group_merged.csv', sep=',')
        
        pdt.assert_frame_equal(df, df_compare)
        
        
    def test_extract_M_and_L(self):
        protein_groups = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_group_merged.csv', sep=',')
        precursor_rollup = DynamicSilac('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/', protein_groups)
        
        
        df = precursor_rollup.extract_M_and_L(protein_groups)
        
        # df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_groups.csv', sep=',', index=False)
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/protein_groups.csv', sep=',')
        
        pdt.assert_frame_equal(df, df_compare)
        
        
if __name__ == '__main__':
    unittest.main()