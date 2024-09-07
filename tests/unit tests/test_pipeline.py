# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 07:43:07 2024

@author: robbi
"""

import unittest
from silac_dia_tools.workflow.pipeline import Pipeline as pipeline


class TestPipeline(unittest.TestCase):
 
    def test_pipeline_initialization(self):
        # Initialize the pipeline with the mocked file paths and parameters
        test_no_spike = pipeline('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/', 'test_params.json', method='dynamic_silac_dia', pulse_channel='M', metadata_file='meta.csv')
        test_spike = pipeline('C:/phd projects/silac_dia_tools_1.0/tests/unit tests/unit test data/spike/', 'test_params.json', method='dynamic_dia_sis', pulse_channel='M', metadata_file='meta.csv')

        # Check if the pipeline object is initialized correctly
        self.assertEqual(test_no_spike.method, 'dynamic_silac_dia')
        self.assertEqual(test_no_spike.pulse_channel, 'M')
        
        self.assertEqual(test_spike.method, 'dynamic_dia_sis')
        self.assertEqual(test_spike.pulse_channel, 'M')
        
   
    

if __name__ == '__main__':
    unittest.main()