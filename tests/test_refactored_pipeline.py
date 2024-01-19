
from icecream import ic
from silac_dia_tools.pipeline.pipeline_refactored import Pipeline as pileline
import pandas as pd 

'''attempting to format the silac channels first then filter afterwards. Filter columns to keep in this step are:
    Parameters used:
        'Lib.PG.Q.Value'
Precursor.Charge > 1
Mass.Evidence > 0.5
Global.PG.Q.Value < 0.01
Channel.Q.Value < 0.03
Translated.Q.Value < 0.03
Translated.Quality >= 0.05

additional columns may be required


'''
if __name__ == "__main__":
    
 
   
    
    path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H refactored/'
    # path = 'G:/My Drive/Data\data/240112 poc4 test/new pipeline and stats/'
    pipeline = pileline( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta='meta.csv')
    # pipeline.make_metadata()
    pipeline.preprocess_href() # in href mode
    pipeline.generate_reports()
    
    # -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:23:34 2024

@author: robbi
"""

