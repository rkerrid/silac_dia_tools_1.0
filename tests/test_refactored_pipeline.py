
from icecream import ic
from silac_dia_tools.pipeline_refactored.pipeline import Pipeline as pileline
import pandas as pd 
pd.set_option('display.max_columns', None)
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

def test_formatting(df):
    
    # Pivot the DataFrame
    pivot_df = df.pivot_table(index=['Run', 'Protein.Group', 'Precursor'], 
                              columns='Label', 
                              values=['Precursor.Quantity', 'Ms1.Translated', 'Precursor.Translated', 'Precursor.Charge', 'Mass.Evidence', 'Global.PG.Q.Value', 'Channel.Q.Value', 'Translated.Q.Value', 'Translated.Quality', 'Lib.PG.Q.Value'],
                              aggfunc='first')
    
    # Flatten the MultiIndex in columns
    pivot_df.columns = [' '.join(col).strip() for col in pivot_df.columns.values]
    
    # Reset index to turn them into columns
    pivot_df.reset_index(inplace=True)
    print('after formatting')
    print(df.head())
    return df

if __name__ == "__main__":
    
 
   
    
    path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H refactored/'
    # path = 'G:/My Drive/Data\data/240112 poc4 test/new pipeline and stats/'
    pipeline = pileline( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta='meta.csv')
    # pipeline.make_metadata()
    df = pipeline.preprocessor.import_report()
    print(df.head())
    test_formatting(df)
    # pipeline.preprocess_pipeline() # in href mode5
    # pipeline.generate_reports()
    
    # -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:23:34 2024

@author: robbi
"""

# try this later 

# # Pivot for each label
# pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor'], aggfunc='first').add_suffix(' L')
# pivot_M = df[df['Label'] == 'M'].pivot_table(index=['Run', 'Protein.Group', 'Precursor'], aggfunc='first').add_suffix(' M')
# pivot_H = df[df['Label'] == 'H'].pivot_table(index=['Run', 'Protein.Group', 'Precursor'], aggfunc='first').add_suffix(' H')

# # Merge the pivoted DataFrames
# merged_df = pd.concat([pivot_L, pivot_M, pivot_H], axis=1)

# # Reset index to make 'Run', 'Protein.Group', and 'Precursor' as columns
# merged_df.reset_index(inplace=True)

# # Display the reformatted DataFrame
# print(merged_df)