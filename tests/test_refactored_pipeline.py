
from icecream import ic
from silac_dia_tools.pipeline_refactored.pipeline import Pipeline as pileline
from silac_dia_tools.pipeline_refactored import generate_protein_groups
# from silac_dia_tools.pipeline_refactored import calculate_intensities_r
# from silac_dia_tools.pipeline_refactored import precursor_rollup  
import pandas as pd 
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np
import operator
import json
import os
import seaborn as sns 


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


if __name__ == "__main__":


    path = 'G:/My Drive/Data/data/20240216 eIF3D timsTOF/'
    # path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/'
    # path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/3d G3 G2/'
    path = 'G:/My Drive/Data/data/20240306 eIF 5 lines/G1 E/'
    # path = 'G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/N/'
    path = 'G:/My Drive/Data/data/240112 poc4 test/20240314 adapted pipeline/H/'
    # path = 'G:/My Drive/Data/main experiments/20240219 baby benchmark for pydia_sis/'
    pipeline = pileline( f'{path}', 'test_params.json', contains_reference = True, method = 'dynamic_dia_sis', pulse_channel="M", meta='meta.csv')
    # pipeline.make_metadata()

    pipeline.execute_pipeline()
 
  
    
    
    
    
    # df, filtered_out, contaminants = pipeline.preprocessor.import_report() 
    # # counts_barplot(df, 'L & M loose filtering')
    
    # # precursor and protein groups
    # # generate precursor report.tsv with precursor.quantity replaced by summing up precursor tranlsated vals
    # df = generate_protein_groups.format_silac_channels(df)
    # # split data (at this point will only work with pre_df and work on median of bot ms1 and pre later)
    # pre_df, ms1_df = generate_protein_groups.split_data_by_intensity_type(df)
    # # calculate precursor ratios and format for input into dlfq
    # pre_df = generate_protein_groups.calculate_precursor_ratios(pre_df, 'Precursor.Translated')
    
    # # calculate protein groups without normalizing by either method
    # protein_groups_unnormalized = generate_protein_groups.compute_protein_level(pre_df)
    
    # # Calculate intensities
    # href_df = calculate_intensities_r.calculate_href_intensities(protein_groups_unnormalized)
    
    # # output href and unnormalized protein groups.csv
    # dfs = calculate_intensities_r.output_protein_groups(href_df, 'href', path)
    
    
    
    
    # ###### DLFQ pipeline POC
    # path = 'G:/My Drive/Data/data/poc4/N/'
    # ### Try with DLFQ
    # pipeline = pileline( f'{path}', 'test_params.json', contains_reference = False, pulse_channel="M", meta='meta.csv')
    # df, filtered_out, contaminants = pipeline.preprocessor.import_report() 
    
    # # precursor and protein groups
    # # generate precursor report.tsv with precursor.quantity replaced by summing up precursor tranlsated vals
    # df = generate_protein_groups.format_silac_channels_dlfq(df)
    # # split data (at this point will only work with pre_df and work on median of bot ms1 and pre later)
    # pre_df, ms1_df = generate_protein_groups.split_data_by_intensity_type_dlfq(df)
    # # calculate precursor ratios and format for input into dlfq
    # pre_df = generate_protein_groups.calculate_precursor_ratios_dlfq(pre_df, 'Precursor.Translated')
    
    # print('writing report tsv')
    # pre_df.to_csv(f'{path}report_.tsv', sep='\t')
    # # calculate protein groups without normalizing by either method
    # silac_precursors_file = f'{path}report_.tsv'
    # dlfq_output_file = f'{path}dlfq_protein_intensities.tsv'
    # from silac_dia_tools.pipeline.utils import dlfq_functions as dlfq

    # # # silac_precursors.to_csv(silac_precursors_file, sep='\t')
    # dlfq.run_lfq(silac_precursors_file, file=dlfq_output_file, num_cores=1)
    
    # protein_groups_unnormalized = generate_protein_groups.compute_protein_level_dlfq(pre_df)
    # dlfq_df = pd.read_csv(dlfq_output_file, sep='\t')
    
    # # Drop the 'Unnamed: 0' column
    # dlfq_df = dlfq_df.drop(columns=['Unnamed: 0', 'protein'])
    
    # # Melt the DataFrame
    # long_df = pd.melt(dlfq_df, id_vars=['Protein.Group'], var_name='Run', value_name='Intensity')
    # merged_dlfq = merge_dlfq_intensities(protein_groups_unnormalized, long_df)
    # m,l = output_protein_groups_dlfq(merged_dlfq, 'dlfq', path)





    
    # barplot_after_all_filtering(df, 'eIF4F Pilot (loose post filtering)')
    
    # ecoli_df = df[df['Protein.Group'].str.contains('ECOLI_')]
    # ecoli_df_rep1 = ecoli_df[ecoli_df['Run'].str.contains('_3')]
    
    # human_df = df[df['Protein.Group'].str.contains('HUMAN_')]
    # human_df_rep1 = human_df[human_df['Run'].str.contains('_3')]
    
    
    # barplot_after_all_filtering_benchmark(ecoli_df_rep1, 'Ecoli After filter by H then loose')
    # barplot_after_all_filtering_benchmark(human_df_rep1, 'Human After filter by H then loose')
# You can then inspect duplicates_df to understand which precursors are causing duplicates
    
    
    
    # # List of columns to check for duplicates
    # cols = ['Run', 'Protein.Group', 'Precursor', 'Label',  'Precursor.Quantity','Ms1.Translated','Precursor.Translated']  + ['Global.Q.Value','RT']
    # # Exclude columns used for grouping from the duplicate check
    # cols_to_check = [col for col in cols if col not in ['Run', 'Protein.Group', 'Precursor', 'Label']]
    
    # # Apply the custom function to each group
    # grouped = df.groupby(['Run', 'Protein.Group', 'Precursor', 'Label'])
    # is_duplicate_series = grouped.transform(lambda x: check_duplicates(x, cols_to_check))

    # # Flatten the result
    # is_duplicate_series = is_duplicate_series.iloc[:, 0]

    # # Add the 'is_duplicate' column to the original DataFrame
    # df['is_duplicate'] = is_duplicate_series

    # # Filter the DataFrame
    # df = df[df['is_duplicate'] == True][['Run', 'Protein.Group', 'Precursor', 'Label'] + cols_to_check]
    # control_df = df[df['Run'] == 'Janice_20231223_AJW_HSdia_H_control_I']
    
    
 

    
    
    
    # print(df.columns.values.tolist())
    # df = test_formatting(df)
    # # print(df.head())
    # strict_parameters = get_strict_params('filtering_parameters_strict.json')
    # print(strict_parameters)
    
    # label= 'M'
    # plot_hl(df, label, 'H vs M')
    # label= 'L'
    # plot_hl(df, label, 'H vs L')
    
    # filter by H channel
    
    # df_filt_h, df_out_h = test_filtering_by_label(df, strict_parameters, 'H')
    
    
   
    
    # # filter other channels with loose filters
    # nan_m_before =count_nans_in_columns(df_filt_h, "M")
    # nan_l_before =count_nans_in_columns(df_filt_h, "L")
    # # add nans for m and l thta dont pass loose
    # loose_params = get_loose_params('filtering_parameters.json')
    # df_filt_h = test_filtering_by_label_add_nan(df_filt_h, loose_params, 'M') # potential bug where it doesnt filter the global.pg.q that was added
    # nan_m = count_nans_in_columns(df_filt_h, "M")
    # df_filt_h= test_filtering_by_label_add_nan(df_filt_h, loose_params, 'L')
    # nan_l =count_nans_in_columns(df_filt_h, "L")
    
    # df_complete = drop_filter_cols(df_filt_h)
    # ic(df_complete)
    # # ic(df.columns.values.tolist())
    # df = df_complete.copy()
    # df = href_intensitie(df) 
    # print('finished renaming')
    # print(df.columns.values.tolist())
    # # ic(df)
    # df = assign_href_intensities(df)
    # cols = ['Run', 'Protein.Group', 'Lib.PG.Q.Value', 'M to stack ratio', 'L to stack ratio', 'H normalized total intensity','Total intensity']
    # df = df[cols]
    # ic(df)
    # create_combined_bar_plot(df)
    
    
    
    # label= 'M'
    # plot_hl(df, label, 'H vs M')
    # label= 'L'
    # plot_hl(df, label, 'H vs L')
    # # roll up to PG
    # protein_groups = calculate_protein_level_ratios(df_filt_all)
    # ic(protein_groups)
    
    
    
    # pipeline.preprocess_pipeline() # in href mode5
    # pipeline.generate_reports()
    
    # -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 16:23:34 2024

@author: robbi
"""

