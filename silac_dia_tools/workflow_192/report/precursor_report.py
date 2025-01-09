
# """
# Created on Sun Sep  8 14:41:41 2024

# @author: robbi
# """
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
import logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)


def create_precursor_report(path):
    input_path = f'{path}preprocessing/'
    output_path = f'{path}/reports/'
    
    
    df = pd.read_csv(f'{input_path}precursors.csv', sep=',')

    # Initialize a PDF file
    with PdfPages(f'{output_path}precursor_report.pdf') as pdf:
        # Call the functions and pass the PDF object
        plot_summed_intensity(df, filtered=True, pdf=pdf)
        plot_id_counts(df, filtered=True, pdf=pdf)
    
    
    
def plot_summed_intensity(df, filtered=False, pdf=None):
    df_filtered = df
    # if filtered:
    #     df_filtered = df[df['filter_passed_H'].astype('bool')]

    metrics = ['precursor_quantity', 'precursor_translated', 'ms1_translated']
    suffixes = ['_L', '_M', '_H', '_pulse']
    suffix_labels = {'_L': 'Light', '_M': 'Medium', '_H': 'Heavy', '_pulse': 'Pulse'}

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))  # Increase the width to accommodate the legend

    for i, metric in enumerate(metrics):
        # Dynamically find which suffixes are available for this metric
        cols = [metric + suffix for suffix in suffixes if metric + suffix in df_filtered.columns]
        
        if not cols:  # Skip if no columns are found for this metric
            continue

        sums = df_filtered.groupby('Run')[cols].sum()
        sums = sums.reset_index()
        sums_melted = sums.melt(id_vars=['Run'], value_vars=cols, var_name='Label', value_name='Sum')
        sums_melted['Suffix'] = sums_melted['Label'].str.extract(r'(_L|_M|_H|_pulse)$')

        unique_suffixes_in_data = sums_melted['Suffix'].dropna().unique()

        # Create a dynamic palette based on the number of unique suffixes in the current plot
        dynamic_palette = [sns.color_palette('Set2', len(unique_suffixes_in_data))[list(unique_suffixes_in_data).index(suffix)] for suffix in unique_suffixes_in_data]

        sns.barplot(x='Run', y='Sum', hue='Suffix', data=sums_melted, palette=dynamic_palette, ax=axes[i])
        axes[i].set_title(f'Sum of Intensities for {metric}', fontsize=14)
        axes[i].set_xlabel('Run', fontsize=12)
        axes[i].tick_params(axis='x', rotation=45)

    # Create the legend based on the suffixes that were present in the dataset
    handles = [mpatches.Patch(color=dynamic_palette[list(unique_suffixes_in_data).index(suffix)], label=suffix_labels[suffix]) for suffix in unique_suffixes_in_data]
    
    # Move the legend closer and adjust the layout
    fig.legend(handles=handles, loc='center right', title='SILAC channel', bbox_to_anchor=(1.05, 0.5), borderaxespad=0)

    plt.tight_layout()

    # Save the plot to the PDF with tight bounding box to include the legend
    if pdf:
        pdf.savefig(fig, bbox_inches='tight')  # 'tight' ensures everything, including the legend, is saved
        plt.close(fig)
    else:
        plt.show()


def plot_id_counts(df, filtered=False, pdf=None):
    df_filtered = df
    # if filtered:
    #     df_filtered = df[df['filter_passed_H'].astype('bool')]

    metrics = ['precursor_quantity', 'precursor_translated', 'ms1_translated']
    suffixes = ['_L', '_M', '_H', '_pulse']
    suffix_labels = {'_L': 'Light', '_M': 'Medium', '_H': 'Heavy', '_pulse': 'Pulse'}

    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)  # Increase the width to accommodate the legend

    for i, metric in enumerate(metrics):
        # Dynamically find which suffixes are available for this metric
        cols = [metric + suffix for suffix in suffixes if metric + suffix in df_filtered.columns]
        
        if not cols:  # Skip if no columns are found for this metric
            continue

        counts = df_filtered.groupby('Run')[cols].apply(lambda x: (x != 0).sum())
        counts = counts.reset_index()
        counts_melted = counts.melt(id_vars=['Run'], value_vars=cols, var_name='Label', value_name='Count')
        counts_melted['Suffix'] = counts_melted['Label'].str.extract(r'(_L|_M|_H|_pulse)$')

        unique_suffixes_in_data = counts_melted['Suffix'].dropna().unique()

        # Create a dynamic palette based on the number of unique suffixes in the current plot
        dynamic_palette = [sns.color_palette('Set2', len(unique_suffixes_in_data))[list(unique_suffixes_in_data).index(suffix)] for suffix in unique_suffixes_in_data]

        sns.barplot(x='Run', y='Count', hue='Suffix', data=counts_melted, palette=dynamic_palette, ax=axes[i])
        axes[i].set_title(f'Non-Zero Counts for {metric}', fontsize=14)
        axes[i].set_xlabel('Run', fontsize=12)
        axes[i].tick_params(axis='x', rotation=45)

    # Create the legend based on the suffixes that were present in the dataset
    handles = [mpatches.Patch(color=dynamic_palette[list(unique_suffixes_in_data).index(suffix)], label=suffix_labels[suffix]) for suffix in unique_suffixes_in_data]
    
    # Move the legend closer and adjust the layout
    fig.legend(handles=handles, loc='center right', title='SILAC channel', bbox_to_anchor=(1.05, 0.5), borderaxespad=0)

    plt.tight_layout()

    # Save the plot to the PDF with tight bounding box to include the legend
    if pdf:
        pdf.savefig(fig, bbox_inches='tight')  # 'tight' ensures everything, including the legend, is saved
        plt.close(fig)
    else:
        plt.show()

if __name__ == '__main__':
    path = r'G:\My Drive\Data\analysis\20240824 BM pilot 2 vs 3 channel analysis\astral\data\no spike\preprocessing\\'
    df = pd.read_csv(f'{path}precursors.csv', sep=',')

    # Initialize a PDF file
    with PdfPages(f'{path}output_plots2.pdf') as pdf:
        # Call the functions and pass the PDF object
        plot_summed_intensity(df, filtered=True, pdf=pdf)
        plot_id_counts(df, filtered=True, pdf=pdf)
