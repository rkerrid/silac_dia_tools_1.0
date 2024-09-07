# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 14:39:00 2024

@author: robbi
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 08:55:46 2024

@author: robbi
"""


import unittest
from silac_dia_tools.workflow.dynamic_dia_sis import DynamicDiaSis
import pandas as pd
import pandas.testing as pdt
from icecream import ic

class TestDynamicDiaSis(unittest.TestCase):
 
    def test_initialize_class(self):
        report = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/imported_reformatted_spike.csv', sep=',')
        precursor_rollup = DynamicDiaSis('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/', report)
    
    def test_get_filtered_set(self):
        report = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/imported_reformatted_spike.csv', sep=',')
        precursor_rollup = DynamicDiaSis('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/', report)
        
        df = precursor_rollup.get_filtered_set(report)
       
        # df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/filtered_spike.csv', index=False)
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/filtered_spike.csv', sep=',')
      
        
        self.assertEqual(len(df), len(df_compare), "DataFrames should have the same number of rows")
        
    def test_generate_href_df(self):
        df = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/filtered_spike.csv', sep=',')
        precursor_rollup = DynamicDiaSis('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/', df)
        
        href = precursor_rollup.generate_href_df(df)
        # href.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/href.csv', index=False)
        
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/href.csv', sep=',')
        
        pdt.assert_frame_equal(href, df_compare)
        
    def test_get_precursor_ratios(self):
        df = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/filtered_spike.csv', sep=',')
        precursor_rollup = DynamicDiaSis('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/', df)
        
        df_l, df_m = precursor_rollup.get_precursor_ratios(df)
        df_l = df_l.reset_index(drop=True)
        df_m = df_m.reset_index(drop=True)
        
        # df_l.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/precursor_ratios_l.csv', index=False)
        # df_m.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/precursor_ratios_m.csv', index=False)
        
        df_compare_l = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/precursor_ratios_l.csv', sep=',')
        df_compare_m = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/precursor_ratios_m.csv', sep=',')
        
        pdt.assert_frame_equal(df_l, df_compare_l)
        pdt.assert_frame_equal(df_m, df_compare_m)

    def test_get_protein_ratios(self):
        df = 0
        precursor_rollup = DynamicDiaSis('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/', df)
        
        df_l = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/precursor_ratios_l.csv', sep=',')
        df_m = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/precursor_ratios_m.csv', sep=',')
        
        df_l = precursor_rollup.get_protein_ratios(df_l, 'L')
        df_m = precursor_rollup.get_protein_ratios(df_m, 'M')
        df_l = df_l.reset_index(drop=True)
        df_m = df_m.reset_index(drop=True)
        
        # df_l.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/protein_ratios_l.csv', index=False)
        # df_m.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/protein_ratios_m.csv', index=False)
        
        df_compare_l = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/protein_ratios_l.csv', sep=',')
        df_compare_m = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/protein_ratios_m.csv', sep=',')
        
        pdt.assert_frame_equal(df_l, df_compare_l)
        pdt.assert_frame_equal(df_m, df_compare_m)
        
    def test_merge_data(self):
        df = 0
        precursor_rollup = DynamicDiaSis('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/', df)
        
        pg_l = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/protein_ratios_l.csv', sep=',')
        pg_m = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/protein_ratios_m.csv', sep=',')
        href = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/href.csv', sep=',')
        
        df = precursor_rollup.merge_data(href, pg_l, pg_m )
        
        # df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/protein_groups.csv', index=False)
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/protein_groups.csv', sep=',')
        
        pdt.assert_frame_equal(df, df_compare)



        
        
    
        
        
if __name__ == '__main__':
    unittest.main()