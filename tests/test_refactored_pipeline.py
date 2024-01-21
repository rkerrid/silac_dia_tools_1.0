
from icecream import ic
from silac_dia_tools.pipeline_refactored.pipeline import Pipeline as pileline
# from silac_dia_tools.pipeline_refactored import precursor_rollup  
import pandas as pd 
import matplotlib.pyplot as plt
import tqdm
import numpy as np
import operator
import json
import os

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
def get_strict_params(parameter_file):
    json_path = f"C:/phd projects/silac_dia_tools1.0/silac_dia_tools/configs/{parameter_file}"
    with open(json_path, 'r') as file:
        return json.load(file)

def get_loose_params(parameter_file):
    json_path = f"C:/phd projects/silac_dia_tools1.0/silac_dia_tools/configs/{parameter_file}"
    with open(json_path, 'r') as file:
        return json.load(file)



def test_formatting(df):
    
   # try this later 
   
   # Pivot for each label
   pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor'], aggfunc='first').add_suffix(' L')
   # ic(pivot_L)
   pivot_M = df[df['Label'] == 'M'].pivot_table(index=['Run', 'Protein.Group', 'Precursor'], aggfunc='first').add_suffix(' M')
   # ic(pivot_M)
   pivot_H = df[df['Label'] == 'H'].pivot_table(index=['Run', 'Protein.Group', 'Precursor'], aggfunc='first').add_suffix(' H')
   # ic(pivot_H)
   # Merge the pivoted DataFrames
   merged_df = pd.concat([pivot_L, pivot_M, pivot_H], axis=1)
   # print('before pivot')
   # ic(df)
   # Reset index to make 'Run', 'Protein.Group', and 'Precursor' as columns
   merged_df.reset_index(inplace=True)

   # Display the reformatted DataFrame
   return merged_df

def test_filtering_by_label(df, params, label):
    # initialize boolen mask the size of the df to be switched to true or false when looping through columns and conditions in params
    filtering_condition = pd.Series([True] * len(df), index=df.index)
    ops = {
        "==": operator.eq, "<": operator.lt, "<=": operator.le,
        ">": operator.gt, ">=": operator.ge
    }

    # Loop through the filtering columns and append ' H'
    for column, condition in params['apply_filters'].items():
        modified_column = column + f' {label}'  # Append ' H' to the column name
        if modified_column in df.columns:  # Check if the modified column exists in the DataFrame
            op = ops[condition['op']]
            #&= means if the row hasnt passed a previous filter it will remain False
            filtering_condition &= op(df[modified_column], condition['value'])

    filtered_chunk = df[filtering_condition]
    chunk_filtered_out = df[~filtering_condition]
    return filtered_chunk, chunk_filtered_out


import numpy as np

def test_filtering_by_label_add_nan(df, params, label):
    filtering_condition = pd.Series([True] * len(df), index=df.index)
    ops = {
        "==": operator.eq, "<": operator.lt, "<=": operator.le,
        ">": operator.gt, ">=": operator.ge
    }

    for column, condition in params['apply_filters'].items():
        modified_column = f'{column} {label}'
        if modified_column in df.columns:
            op = ops[condition['op']]
            filtering_condition &= op(df[modified_column], condition['value'])
            filtered_chunk = df.copy()
            nan_cols = [f'Precursor.Translated {label}', f'Precursor.Quantity {label}', f'Ms1.Translated {label}']
            for col in nan_cols:
                if col in df.columns:
                    filtered_chunk.loc[~filtering_condition, col] = np.nan

    return filtered_chunk

def count_nans_in_columns(df, labels):
    # Columns to check for NaNs, adjust these as per your DataFrame's structure
    columns_to_check = [f'Precursor.Translated {label}' for label in labels] + \
                       [f'Precursor.Quantity {label}' for label in labels] + \
                       [f'Ms1.Translated {label}' for label in labels]

    nan_counts = {}
    for col in columns_to_check:
        if col in df.columns:
            nan_counts[col] = df[col].isna().sum()
    
    return nan_counts


# def plot_data(df, label, title):
#     # Plotting scatter plot
#     plt.figure(figsize=(10, 6))
#     plt.scatter(df[f'Translated.Quality {label}'], df[f'Precursor.Translated {label}'], color='blue')
#     plt.scatter(df[f'Translated.Q.Value {label}'], df[f'Precursor.Translated {label}'], color='red')
#     plt.xlabel(f'Translated.Quality {label}')
#     plt.ylabel(f'Precursor.Translated {label}')
#     plt.title(f'Scatter Plot of Precursor.Translated {label} vs Translated.Quality {label} {title}')
#     plt.grid(True)
#     plt.legend()
#     plt.show()

def plot_translated_qval(df, label, title):
    # Plotting scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df[f'Translated.Q.Value {label}'], df[f'Precursor.Translated {label}'], color='red', label=f'Translated.Q.Value {label}')
    plt.xlabel(f'Translated.Q val {label}')
    plt.ylabel(f'Precursor.Translated {label}')
    plt.title(f'Scatter Plot of Precursor.Translated {label} vs Translated.Q val {label} {title}')
    plt.grid(True)
    plt.legend()  # This will display the legend with correct labeling of colors
    plt.show()

def plot_translated_quality(df, label, title):
    # Plotting scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df[f'Translated.Quality {label}'], df[f'Precursor.Translated {label}'], color='blue', label=f'Translated.Quality {label}')
    plt.xlabel(f'Translated.Quality {label}')
    plt.ylabel(f'Precursor.Translated {label}')
    plt.title(f'Scatter Plot of Precursor.Translated {label} vs Translated.Quality {label} {title}')
    plt.grid(True)
    plt.legend()  # This will display the legend with correct labeling of colors
    plt.show()



# def apply_filter_conditions(df, params, suffix):
#     filtering_condition = pd.Series([True] * len(df), index=df.index)
#     ops = {
#         "==": operator.eq, "<": operator.lt, "<=": operator.le,
#         ">": operator.gt, ">=": operator.ge
#     }
#     for column, condition in params['apply_filters'].items():
#         modified_column = column + suffix
#         if modified_column in df.columns:
#             op = ops[condition['op']]
#             filtering_condition &= op(df[modified_column], condition['value'])
#     return filtering_condition

# def test_filtering_formatted(df, loose_params):
#     print("called loose params filtering post heavy ")

#     filtering_condition_M = apply_filter_conditions(df, loose_params, 'M')
#     filtering_condition_L = apply_filter_conditions(df, loose_params, 'L')

#     # Combine the conditions
#     #must both be true to keep row
#     combined_condition = filtering_condition_M & filtering_condition_L

#     filtered_chunk = df[combined_condition]
#     chunk_filtered_out = df[~combined_condition]

#     return filtered_chunk, chunk_filtered_out



def _select_ms1_translated(group):
    return group[group['quantity type'] == 'Ms1.Translated']

def calculate_protein_level_ratios(df):
    print("Calculating ratios from precursor information")
    protein_precursors = df.groupby(['Run', 'Protein.Group'])
    ic(protein_precursors)
    protein_data = []
    protein_count, protein_missed = 0, 0

    for name, group in protein_precursors:
        
        if len(protein_precursors) > 2:
            # median_ratios = _calculate_median_log_ratios(ms1_group)
            total_intensity = np.sum(df['Precursor.Quantity L'])
            new_row = _create_protein_row(group, total_intensity)
            protein_data.append(new_row)
            protein_count += 1
        else:
            protein_missed += 1

    _print_protein_counts(protein_count, protein_missed, len(protein_precursors))
    return pd.DataFrame(protein_data)

def _calculate_median_log_ratios( group):
    median_log2_ratios = np.median(np.log2(group[['L to stack ratio', 'M to stack ratio', 'H to stack ratio']]), axis=0)
    return np.exp2(median_log2_ratios)

def _create_protein_row(group, total_intensity):
    return {
        'Run': group['Run'].iloc[0],
        'Protein.Group': group['Protein.Group'].iloc[0],
        'Total intensity': total_intensity,
        # 'L to stack ratio': median_ratios[0],
        # 'M to stack ratio': median_ratios[1],
        # 'H to stack ratio': median_ratios[2],
        # 'L intensity': median_ratios[0] * total_intensity,
        # 'M intensity': median_ratios[1] * total_intensity,
        # 'H intensity': median_ratios[2] * total_intensity
    }

def _print_protein_counts(protein_count, protein_missed, total):
    print(f'Total proteins counted: {protein_count}')
    print(f'Total sets of precursors that didn\'t meet minimum unique precursor requirements: {protein_missed} out of {total}')
    
    
def plot_hl(df, label, title):
    # Plotting scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df[f'Precursor.Translated H'],df[f'Precursor.Translated {label}'], color='blue', label=f'precursor translated {label}')
    plt.ylabel(f'Precursor.Translated {label}')
    plt.xlabel(f'Precursor.Translated H')
    plt.title(f'Scatter Plot of Precursor.Translated H vs pre translated {label} {title}')
    # plt.xlim(1000000000)
    # plt.ylim(1000000000)
    plt.grid(True)
    plt.legend()  # This will display the legend with correct labeling of colors
    plt.show()
    
def drop_filter_cols(df):

    cols =  ['Run','Protein.Group','Ms1.Translated L','Ms1.Translated M','Ms1.Translated H','Lib.PG.Q.Value H']
    df = df[cols]

    old_names = ['Ms1.Translated L', 'Ms1.Translated M', 'Ms1.Translated H', 'Lib.PG.Q.Value H']
    new_names = ['L intensity', 'M intensity', 'H intensity', 'Lib.PG.Q.Value']
    
    # Create a mapping from old names to new names
    name_mapping = dict(zip(old_names, new_names))
    
    # Rename the columns
    df_copy = df.copy()
    df_copy.rename(columns=name_mapping, inplace=True)
    # print('renamed')
    return df_copy

def calculate_ratios(df):
    df['Precursor.Quantity'] = df['H intensity'] + df['M intensity'] + df['L intensity']
    for label in ['L', 'M', 'H']:
        ratio_column = f'{label} to stack ratio'
        df[ratio_column] = df[f'{label} intensity'] / df['Precursor.Quantity']
        # ic(df)
        print("after calculating ratios")
    return df


def href_intensitie(df):
    print('Begin href intensities')

    # Calculate the median 'H intensity' and rename it to 'h_ref'
    median_h_ref = df.groupby('Protein.Group')['H intensity'].median().reset_index().rename(columns={'H intensity': 'h_ref'})
    df = pd.merge(df, median_h_ref, on='Protein.Group', how='left')
    df['H to stack ratio'] = df['H intensity'] / df['h_ref']
    df['M to stack ratio'] = df['M intensity'] / df['h_ref']
    df['L to stack ratio'] = df['L intensity'] / df['h_ref']
    print('after merge')
    print('Complete href intensities')
    return df

def assign_href_intensities(df):
    df['H normalized total intensity'] = df['h_ref'] / df['H to stack ratio']
    df['M'] = df['H normalized total intensity'] *  df['M to stack ratio']
    df['L'] = df['H normalized total intensity']  *  df['L to stack ratio']
    # ic(df)
    for channel in ['H', 'M', 'L']:
        df[f'{channel} intensity'] = df['H normalized total intensity'] * df[f'{channel} to stack ratio']
    print('after asugning differbte ubtebsutues after birnakuzatuib')
    # ic(df)
    df['Total intensity'] = df['L intensity'] + df['M intensity']
    # df['NSP intensity'] = df['M intensity']
            
    return df
    

def create_combined_bar_plot(df):
    columns = ['M to stack ratio', 'L to stack ratio', 'H normalized total', 'Total intensity']
    
    # Calculate the total of non-NaN values for each column
    totals = [df[col].dropna().sum() if col in df.columns else 0 for col in columns]

    # Create a bar plot for these totals
    plt.figure(figsize=(10, 6))
    plt.bar(columns, totals, color=['blue', 'green', 'red', 'orange'])
    plt.title('Total of Non-NaN Values for Each Column')
    plt.ylabel('Total Value')
    plt.xlabel('Column')

    plt.show()

if __name__ == "__main__":
    
 
   
    
    path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H refactored/'
    # path = 'G:/My Drive/Data\data/240112 poc4 test/new pipeline and stats/'
    pipeline = pileline( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta='meta.csv')
    # pipeline.make_metadata()
    df = pipeline.preprocessor.import_report()
    # print(df.columns.values.tolist())
    df = test_formatting(df)
    # print(df.head())
    strict_parameters = get_strict_params('filtering_parameters_strict.json')
    # print(strict_parameters)
    
    # label= 'M'
    # plot_hl(df, label, 'H vs M')
    # label= 'L'
    # plot_hl(df, label, 'H vs L')
    
    # filter by H channel
    
    df_filt_h, df_out_h = test_filtering_by_label(df, strict_parameters, 'H')
    
    
    # label= 'M'
    # plot_hl(df, label, 'H vs M')
    # label= 'L'
    # plot_hl(df, label, 'H vs L')
    # ic(df_filt)
    # ic(df_out)
    # plot_translated_qval(df,'H', 'before df')
    # plot_translated_qval(df_filt, 'H', 'filtered df')
    # plot_translated_qval(df_out,'H', 'filtered out df')
    
    # plot_translated_qval(df,'L', 'before df')
    # plot_translated_qval(df_filt, 'L', 'filtered df')
    # plot_translated_qval(df_out,'L', 'filtered out df')
    
    # plot_translated_qval(df,'M', 'before df')
    # plot_translated_qval(df_filt, 'M', 'filtered df')
    # plot_translated_qval(df_out,'M', 'filtered out df')
    
    
    
    # plot_translated_quality(df,'H', 'before df')
    # plot_translated_quality(df_filt, 'H', 'filtered df')
    # plot_translated_quality(df_out,'H', 'filtered out df')
    
    # plot_translated_quality(df,'L', 'before df')
    # plot_translated_quality(df_filt, 'L', 'filtered df')
    # plot_translated_quality(df_out,'L', 'filtered out df')
    
    # plot_translated_quality(df,'M', 'before df')
    # plot_translated_quality(df_filt, 'M', 'filtered df')
    # plot_translated_quality(df_out,'M', 'filtered out df')
    
    # filter other channels with loose filters
    nan_m_before =count_nans_in_columns(df_filt_h, "M")
    nan_l_before =count_nans_in_columns(df_filt_h, "L")
    # add nans for m and l thta dont pass loose
    loose_params = get_loose_params('filtering_parameters.json')
    df_filt_h = test_filtering_by_label_add_nan(df_filt_h, loose_params, 'M') # potential bug where it doesnt filter the global.pg.q that was added
    nan_m = count_nans_in_columns(df_filt_h, "M")
    df_filt_h= test_filtering_by_label_add_nan(df_filt_h, loose_params, 'L')
    nan_l =count_nans_in_columns(df_filt_h, "L")
    
    df_complete = drop_filter_cols(df_filt_h)
    ic(df_complete)
    # ic(df.columns.values.tolist())
    df = df_complete.copy()
    df = href_intensitie(df) 
    print('finished renaming')
    print(df.columns.values.tolist())
    # ic(df)
    df = assign_href_intensities(df)
    cols = ['Run', 'Protein.Group', 'Lib.PG.Q.Value', 'M to stack ratio', 'L to stack ratio', 'H normalized total intensity','Total intensity']
    df = df[cols]
    ic(df)
    create_combined_bar_plot(df)
    
    # print(df.columns.values.tolist())
    # plot_translated_qval(df_filt,'H', 'before (after h filt loose on other channels) df')
    # plot_translated_qval(df_filt_all, 'H', 'after (after h filt loose on other channels) df')
    # plot_translated_qval(df_filt_out,'H', 'filtered out after (after h filt loose on other channels)')
    
    # plot_translated_qval(df_filt,'L', 'before (after h filt loose on other channels) df')
    # plot_translated_qval(df_filt_all, 'L', 'after (after h filt loose on other channels) df')
    # plot_translated_qval(df_filt_out,'L', 'filtered out after (after h filt loose on other channels)')
    
    # plot_translated_qval(df_filt,'M', 'before (after h filt loose on other channels) df')
    # plot_translated_qval(df_filt_all, 'M', 'after (after h filt loose on other channels) df')
    # plot_translated_qval(df_filt_out,'M', 'filtered out after (after h filt loose on other channels)')
   
    
    
    # plot_translated_quality(df_filt,'H', 'before (after h filt loose on other channels) df')
    # plot_translated_quality(df_filt_all, 'H', 'after (after h filt loose on other channels) df')
    # plot_translated_quality(df_filt_out,'H', 'filtered out after (after h filt loose on other channels)')
    
    # plot_translated_quality(df_filt,'L', 'before (after h filt loose on other channels) df')
    # plot_translated_quality(df_filt_all, 'L', 'after (after h filt loose on other channels) df')
    # plot_translated_qval(df_filt_out,'L', 'filtered out after (after h filt loose on other channels)')
    
    # plot_translated_quality(df_filt,'M', 'before (after h filt loose on other channels) df')
    # plot_translated_quality(df_filt_all, 'M', 'after (after h filt loose on other channels) df')
    # plot_translated_quality(df_filt_out,'M', 'filtered out after (after h filt loose on other channels)') 
    
    
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

