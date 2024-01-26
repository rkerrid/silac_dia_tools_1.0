
from icecream import ic
from silac_dia_tools.pipeline_refactored.pipeline import Pipeline as pileline
from silac_dia_tools.pipeline_refactored import generate_protein_groups
from silac_dia_tools.pipeline_refactored import calculate_intensities_r
# from silac_dia_tools.pipeline_refactored import precursor_rollup  
import pandas as pd 
import matplotlib.pyplot as plt
import tqdm
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


# def test_formatting(df):
    
#    # try this later 
   
#    # Pivot for each label
#    pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' L')
#    pivot_M = df[df['Label'] == 'M'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' M')
#    pivot_H = df[df['Label'] == 'H'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' H')
#    # Merge the pivoted DataFrames
#    merged_df = pd.concat([pivot_L, pivot_M, pivot_H], axis=1)
#    # Reset index to make 'Run', 'Protein.Group', and 'Precursor' as columns
#    merged_df.reset_index(inplace=True)

#    # remove all rows that contain a NaN under the Label H column (i.e. no H precurosr is present for that row)
#    if 'Label H' in df.columns:
#        df = merged_df.dropna(subset=['Label H'])
#    else:
#        print("Column 'Label H' not found in DataFrame")
#        print(merged_df.columns.values.tolist())
#        df = merged_df
#    return df

def test_formatting(df):
    # Pivot for each label
    pivot_L = df[df['Label'] == 'L'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' L')
    pivot_M = df[df['Label'] == 'M'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' M')
    pivot_H = df[df['Label'] == 'H'].pivot_table(index=['Run', 'Protein.Group', 'Precursor.Id'], aggfunc='first').add_suffix(' H')
    
    # Merge the pivoted DataFrames
    merged_df = pd.concat([pivot_L, pivot_M, pivot_H], axis=1)
    # Reset index to make 'Run', 'Protein.Group', and 'Precursor.Id' as columns
    merged_df.reset_index(inplace=True)

    # remove all rows that contain a NaN under the Label H column (i.e., no H precursor is present for that row)
    if 'Label H' in merged_df.columns:
        # Apply dropna on merged_df instead of df
        merged_df = merged_df.dropna(subset=['Precursor.Quantity H'])
    else:
        print("Column 'Label H' not found in DataFrame")
        print(merged_df.columns.values.tolist())
    
    return merged_df




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
    

def counts_barplot(df, title):
    
    # Filter out rows where 'Precursor.Translated' is NaN or 0
    filtered_df = df[(df['Precursor.Translated'].notna()) & (df['Precursor.Translated'] != 0)]
    
    # Group by 'Run' and 'Label', and count the non-zero, non-NaN 'Precursor.Translated' values
    grouped_df = filtered_df.groupby(['Run', 'Label'])['Precursor.Translated'].count().reset_index()
    
    # Plot
    plt.figure(figsize=(12, 6))
    barplot = sns.barplot(data=grouped_df, x='Run', y='Precursor.Translated', hue='Label')
    plt.xticks(rotation=90)  # Rotate labels for readability
    plt.xticks(rotation=90)  # Rotate labels for readability
    
    # Iterate over the bars
    for bar in barplot.patches:
        # Using the bar's height to place the label
        barplot.text(bar.get_x() + bar.get_width() / 2., bar.get_height(),
                     int(bar.get_height()),  # The label
                     ha='center', va='bottom',
                     rotation=45)  # Rotate label
    
    plt.title(f'Barplot of H M and L precursor counts: {title} settings')
    plt.ylabel('Count of valid vals (Precursor.Translated)')
    plt.xlabel('Run')
    plt.tight_layout()
    
    save_path = f'G:/My Drive/Data/figures/{title}_Precursor.Translated.jpeg'
    plt.savefig(save_path, format='jpeg', dpi=300)
    
    plt.show()
    
    # Filter out rows where 'Precursor.Translated' is NaN or 0
    filtered_df = df[(df['Ms1.Translated'].notna()) & (df['Ms1.Translated'] != 0)]
    
    # Group by 'Run' and 'Label', and count the non-zero, non-NaN 'Precursor.Translated' values
    grouped_df = filtered_df.groupby(['Run', 'Label'])['Ms1.Translated'].count().reset_index()
    
    # Plot
    plt.figure(figsize=(12, 6))
    barplot = sns.barplot(data=grouped_df, x='Run', y='Ms1.Translated', hue='Label')
    plt.xticks(rotation=90)  # Rotate labels for readability
    plt.xticks(rotation=90)  # Rotate labels for readability
    # Iterate over the bars
    for bar in barplot.patches:
        # Using the bar's height to place the label
        barplot.text(bar.get_x() + bar.get_width() / 2., bar.get_height(),
                     int(bar.get_height()),  # The label
                     ha='center', va='bottom',
                     rotation=45)  # Rotate label
    
    plt.title(f'Barplot of H M and L precursor counts: {title} settings')
    plt.ylabel('Count of valid vals (Ms1.Translated)')
    plt.xlabel('Run')
    plt.tight_layout()
    
    save_path = f'G:/My Drive/Data/figures/{title}_Ms1.Translated.jpeg'
    plt.savefig(save_path, format='jpeg', dpi=300)
    
    plt.show()
    
    # Filter out rows where 'Precursor.Translated' is NaN or 0
    filtered_df = df[(df['Precursor.Quantity'].notna()) & (df['Precursor.Quantity'] != 0)]
    
    # Group by 'Run' and 'Label', and count the non-zero, non-NaN 'Precursor.Translated' values
    grouped_df = filtered_df.groupby(['Run', 'Label'])['Precursor.Quantity'].count().reset_index()
    
    # Plot
    plt.figure(figsize=(12, 6))
    barplot = sns.barplot(data=grouped_df, x='Run', y='Precursor.Quantity', hue='Label')
    plt.xticks(rotation=90)  # Rotate labels for readability
    # Iterate over the bars
    for bar in barplot.patches:
        # Using the bar's height to place the label
        barplot.text(bar.get_x() + bar.get_width() / 2., bar.get_height(),
                     int(bar.get_height()),  # The label
                     ha='center', va='bottom',
                     rotation=45)  # Rotate label
    
    plt.title(f'Barplot of H M and L precursor counts: {title} settings')
    plt.ylabel('Count of valid vals (Precursor.Quantity)')
    plt.xlabel('Run')
    plt.tight_layout()
    
    save_path = f'G:/My Drive/Data/figures/{title}_Precursor.Quantity.jpeg'
    plt.savefig(save_path, format='jpeg', dpi=300)
    
    
    plt.show()

# def barplot_after_all_filtering(df):
    
#     # Calculate non-NaN counts for each 'Run'
#     grouped = df.groupby('Run')[['Precursor.Quantity H', 'Precursor.Quantity M', 'Precursor.Quantity L']].count().reset_index()
    
#     # Plotting
#     fig, ax = plt.subplots(figsize=(10, 6))
    
#     # Locations of the groups on the x-axis
#     x = range(len(grouped))
    
#     # Plot each of the quantities as a separate bar
#     ax.bar(x, grouped['Precursor.Quantity H'], width=0.2, label='H', align='center')
#     ax.bar(x, grouped['Precursor.Quantity M'], width=0.2, label='M', align='edge')
#     ax.bar(x, grouped['Precursor.Quantity L'], width=0.2, label='L', align='edge')
    
#     # Set the labels and titles
#     ax.set_ylabel('Non-NaN Counts')
#     ax.set_title('Non-NaN Counts of Precursor Quantities by Run')
#     ax.set_xticks(x)
#     ax.set_xticklabels(grouped['Run'], rotation=45)
#     ax.legend()
    
#     # Show the plot
#     plt.show()


def barplot_after_all_filtering(df, title):
    # Calculate non-NaN counts for each 'Run' and 'Label'
    counts_h = df.groupby('Run')['Precursor.Quantity H'].count().reset_index(name='Count')
    counts_h['Label'] = 'H'
    
    counts_l = df.groupby('Run')['Precursor.Quantity L'].count().reset_index(name='Count')
    counts_l['Label'] = 'L'
    
    counts_m = df.groupby('Run')['Precursor.Quantity M'].count().reset_index(name='Count')
    counts_m['Label'] = 'M'
    # Combine the counts into a single DataFrame
    combined_counts = pd.concat([counts_h, counts_l, counts_m ]) #, counts_m

    # Plot using Seaborn
    plt.figure(figsize=(12, 6))
    barplot = sns.barplot(data=combined_counts, x='Run', y='Count', hue='Label')

    # Rotate labels for readability
    plt.xticks(rotation=90)

    # Iterate over the bars for labels
    for bar in barplot.patches:
        # Using the bar's height to place the label
        barplot.text(bar.get_x() + bar.get_width() / 2., bar.get_height(),
                     int(bar.get_height()),  # The label
                     ha='center', va='bottom',
                     rotation=45)  # Rotate label

    plt.title(f'Barplot of H, M, and L Precursor Counts: {title}')
    plt.ylabel('Count of Non-NaN Values')
    plt.xlabel('Run')
    plt.legend(title='Label', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    
    
    save_path = f'G:/My Drive/Data/figures/{title}_Precursor.Quantity.jpeg'
    plt.savefig(save_path, format='jpeg', dpi=300)
    
    plt.show()
    
def barplot_after_all_filtering_benchmark(df, title):
    

    # Calculate non-NaN counts for each 'Run' and 'Label'
    counts_h = df.groupby('Run')['Precursor.Quantity H'].count().reset_index(name='Count')
    counts_h['Label'] = 'H'
    
    counts_l = df.groupby('Run')['Precursor.Quantity L'].count().reset_index(name='Count')
    counts_l['Label'] = 'L'
    
    # counts_m = df.groupby('Run')['Precursor.Quantity M'].count().reset_index(name='Count')
    # counts_m['Label'] = 'M'
    # Combine the counts into a single DataFrame
    combined_counts = pd.concat([counts_h, counts_l ]) 

    # Plot using Seaborn
    plt.figure(figsize=(12, 6))
    barplot = sns.barplot(data=combined_counts, x='Run', y='Count', hue='Label')

    # Rotate labels for readability
    plt.xticks(rotation=90)

    # Iterate over the bars for labels
    for bar in barplot.patches:
        # Using the bar's height to place the label
        barplot.text(bar.get_x() + bar.get_width() / 2., bar.get_height(),
                     int(bar.get_height()),  # The label
                     ha='center', va='bottom',
                     rotation=45)  # Rotate label

    plt.title(f'Barplot of H, M, and L Precursor Counts: {title}')
    plt.ylabel('Count of Non-NaN Values')
    plt.xlabel('Run')
    plt.legend(title='Label', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    
    
    save_path = f'G:/My Drive/Data/figures/{title}_Precursor.Quantity.jpeg'
    plt.savefig(save_path, format='jpeg', dpi=300)
    
    plt.show()
        


if __name__ == "__main__":

 
   
    path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/H refactored/'
    path = 'G:/My Drive/Data/data/20240125 bm filter h then loose filter l/'
    path = 'G:/My Drive/Data/data/eIF4F optimization/'
    # path = 'G:/My Drive/Data\data/240112 poc4 test/new pipeline and stats/'
    pipeline = pileline( f'{path}', 'test_params.json', contains_reference = True, pulse_channel="M", meta='meta.csv')
    df, filtered_out, contaminants = pipeline.preprocessor.import_report() 
    # counts_barplot(df, 'L & M loose filtering')
    
    # precursor and protein groups
    # generate precursor report.tsv with precursor.quantity replaced by summing up precursor tranlsated vals
    df = generate_protein_groups.format_silac_channels(df)
    # split data (at this point will only work with pre_df and work on median of bot ms1 and pre later)
    pre_df, ms1_df = generate_protein_groups.split_data_by_intensity_type(df)
    # calculate precursor ratios and format for input into dlfq
    pre_df = generate_protein_groups.calculate_precursor_ratios(pre_df, 'Precursor.Translated')
    # calculate protein groups without normalizing by either method
    protein_groups_unnormalized = generate_protein_groups.compute_protein_level(pre_df)
    
    # Calculate intensities
    href_df = calculate_intensities_r.calculate_href_intensities(protein_groups_unnormalized)
    
    # output href and unnormalized protein groups.csv
    dfs = calculate_intensities_r.output_protein_groups(href_df, 'href', path)
    
    
    
    
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

