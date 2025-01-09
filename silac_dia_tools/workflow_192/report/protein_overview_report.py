# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 10:15:43 2024

@author: rkerrid
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import os
from silac_dia_tools.workflow.utils import manage_directories


def import_data(path):
    path = f'{path}/protein_groups'
    light = pd.read_csv(f'{path}/light.csv', sep=',')
    nsp = pd.read_csv(f'{path}/nsp.csv', sep=',')
    #href = pd.read_csv(f'{path}/href.csv', sep=',')
    return light, nsp

def summarize_data(df):
    df_counts = df.notna().sum()[1:]
    
    df_cols = df.columns.values.tolist()
    df_cols = df_cols[1:]
    
    df_intensity = df.iloc[:,1:].sum(skipna=True)    
    df_intensity = df_intensity.reset_index()
    df_intensity.columns = ['Experiment', 'Total Intensity']
    df_intensity = df_intensity['Total Intensity'].values.tolist()
    
    summarize_df = {'experiment':df_cols,
                     'counts':df_counts,
                      'total_intensity': df_intensity}

    summarize_df = pd.DataFrame(summarize_df)
    return summarize_df

def create_report(path):
    
    # import NSP and M df
    light, nsp = import_data(path)
    
    # summarize data
    light_summarize = summarize_data(light)
    light_summarize['channel'] = 'light'
    nsp_summarize = summarize_data(nsp)
    nsp_summarize['channel'] = 'nsp'
    
    df = pd.concat([light_summarize, nsp_summarize])
    
    pdf_path = f'{path}/reports/protein_group_overview.pdf'
    
    with PdfPages(pdf_path) as pdf:
        # Title and introduction
        plt.figure(figsize=(11, 8))
        plt.axis('off')

        plt.text(0.5, 0.98, "Protein Groups QC Report", ha='center', va='top', fontsize=15, fontweight='bold')
        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()
        
        plt.figure(figsize=(15, 8))
        sns.barplot(data=df, x='experiment', y='counts', hue='channel')
        plt.title('Number of Protein IDs in Each Experiment')
        plt.xlabel('Experiments')
        plt.ylabel('Number of Protein IDs')
        plt.legend(title='SILAC Channels')
        plt.xticks(rotation=-90)
        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()
        
        plt.figure(figsize=(15, 8))
        sns.barplot(data=df, x='experiment', y='total_intensity', hue='channel')
        plt.title('Summed Protein Intensities in Each Experiment')
        plt.xlabel('Experiments')
        plt.ylabel('Summed Intensity')
        plt.legend(title='SILAC Channels')
        plt.xticks(rotation=-90)
        
        pdf.savefig()  # Saves the current figure into the PDF
        plt.close()

if __name__ == '__main__':
    path = 'G:/My Drive/Data/data/20240530 modifying triple silac/'
    create_report(path)