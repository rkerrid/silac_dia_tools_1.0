import seaborn as sns
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import json

from silac_dia_tools.pipeline.utils import manage_directories


""" 
Currently in progress will focus after pg report and read out things like
1. number of IDs
2. plot histogram of values
3. plot histogram of ratios for each
4. any other QC information
"""
def create_report(path, params, method, pulse_channel):
    # Construct the description string
    params_str = "\n".join([f"{key} {item['op']} {item['value']}" for key, item in params['apply_strict_filters'].items()])
    description = f"Parameters used:\n{params_str}"
    
    # read in data
    light_precursors = pd.read_csv(f'{path}/preprocessing/light_precursors.tsv')
    medium_precursors = pd.read_csv(f'{path}/preprocessing/medium_precursors.tsv')

    counts_df = light_precursors['Run'].value_counts()
    
    # Set up the PDF
    # create_reports_directory(path)
    manage_directories.create_directory(path,'reports')
    output_dir = path + '/reports'
    pdf_path = os.path.join(output_dir, 'silac_precursors_report.pdf')
    
    df_grouped = df.groupby('Run')
    
    
    
    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(8, 11))
        plt.axis('off')
        plt.text(0.5, 0.98, "Silac Precursors QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        plt.text(0.5, 0.85, description, ha='center', va='center', wrap=True)

        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()
        

                

        # Display table of counts
        table_data = [["Run", "DF Count"]]
        for run in counts_df.keys():
            table_data.append([run, counts_df.get(run, 0)])
        
        plt.table(cellText=table_data, cellLoc = 'center', loc='center', colWidths=[0.2,0.2])
        plt.axis('off')
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()

        # Get the list of unique runs
        runs = df_grouped.groups.keys()

        # For each run, plot the histograms and save to PDF
        if method == 'dynamic_dia_sis':
            for run in runs:
                plot_histograms_for_run(run, df_grouped, ['Precursor.Translated H', 'Precursor.Translated M', 'Precursor.Translated L'])
                
                pdf.savefig()  # Saves the current figure into the PDF
                plt.close()
        elif method == 'dia_sis':
            for run in runs:
                plot_histograms_for_run(run, df_grouped, ['Precursor.Translated H', 'Precursor.Translated L'])
                
                # )
                pdf.savefig()  # Saves the current figure into the PDF
                plt.close()
        else:
            for run in runs:
                plot_histograms_for_run(run, df_grouped, [f'Precursor.Translated {pulse_channel}', 'Precursor.Translated L'])
                
                # )
                pdf.savefig()  # Saves the current figure into the PDF
                plt.close()
            
             
def plot_histograms_for_run(run, df_grouped, labels):
    """Plot histograms for a given run from multiple dataframes."""
    colors = sns.color_palette("husl", 3)  # Get a color palette

    # Extract the DataFrame for the specific run
    df_run = df_grouped.get_group(run)

    for label, color in zip(labels, colors):
        # Calculate log2 values
        data = np.log2(df_run[label].dropna())  # Add 1 to avoid log2(0)
        
        # Filter out -inf and inf values
        data = data[np.isfinite(data)]
        plt.hist(data, alpha=0.5, label=label, bins=300, color=color)

    plt.title(f'Histogram for {run}')
    plt.xlabel("log2 Intensity")
    plt.ylabel('Frequency')
    plt.legend()

