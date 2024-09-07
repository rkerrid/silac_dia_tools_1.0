# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 08:07:47 2024

@author: robbi
"""


import unittest
from silac_dia_tools.workflow.preprocessor import Preprocessor 
import pandas as pd
import pandas.testing as pdt


class TestPipeline(unittest.TestCase):
 
    def initialize_preprocessors(self):
        no_spike_meta = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/meta.csv', sep=',')
        no_spike_preprocessor = Preprocessor('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/',  'dynamic_silac_dia', 'M', no_spike_meta)
        
        spike_meta = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/meta.csv', sep=',')
        spike_preprocessor = Preprocessor('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/',  'dynamic_dia_sis', 'M', spike_meta)
   
        return spike_preprocessor, no_spike_preprocessor
    
    def test_no_spike_import(self):
        no_spike_meta = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/meta.csv', sep=',')
        no_spike_preprocessor = Preprocessor('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/',  'dynamic_silac_dia', 'M', no_spike_meta)
        
        df, contams = no_spike_preprocessor.import_report()
        
        # df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_no_spike.csv', sep=',', index=False) 
        # contams.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_no_spike_contams.csv', sep=',', index=False) 
        
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_no_spike.csv', sep=',')
        contams_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_no_spike_contams.csv', sep=',') 
        
        pdt.assert_frame_equal(df, df_compare)
        
        pdt.assert_frame_equal(contams, contams_compare)
        
    def test_spike_import(self):
        spike_meta = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/meta.csv', sep=',')
        spike_preprocessor = Preprocessor('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/',  'dynamic_silac_dia', 'M', spike_meta)
        
        df, contams = spike_preprocessor.import_report()
        
        df.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/imported_spike.csv', sep=',', index=False) 
        contams.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/imported_spike_contams.csv', sep=',', index=False) 
        
        df_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/imported_spike.csv', sep=',')
        contams_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/imported_spike_contams.csv', sep=',') 
        
        pdt.assert_frame_equal(df, df_compare)
        
        pdt.assert_frame_equal(contams, contams_compare)
        
    def test_no_spike_reformat_table(self):
        
        no_spike_meta = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/meta.csv', sep=',')
        no_spike_preprocessor = Preprocessor('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/',  'dynamic_silac_dia', 'M', no_spike_meta)
        
        df_imported = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_no_spike.csv', sep=',')
        
        df_reformated = no_spike_preprocessor.reformat_table(df_imported, 'dynamic_silac_dia', 'M')
        
        # df_reformated.to_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_reformatted_no_spike.csv', sep=',', index=False) 
        df_reformated_compare = pd.read_csv('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/no spike/imported_reformatted_no_spike.csv', sep=',')

        pdt.assert_frame_equal(df_reformated, df_reformated_compare)

if __name__ == '__main__':
    unittest.main()