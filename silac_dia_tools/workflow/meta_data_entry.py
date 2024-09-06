# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 13:52:29 2024

@author: rkerrid
"""

import tkinter as tk
from tkinter import filedialog
from pandastable import Table
import dask.dataframe as dd
import os
import pandas as pd


class MetaDataEntry(tk.Frame):
    def __init__(self, path, parent=None):
        super().__init__(parent)
        self.path = path
        self.initialize_ui()
        
    def initialize_ui(self):
        self.main = self.master
        self.main.geometry('600x400+200+100')
        self.main.title('Data Entry Table')
        self.df = self.create_table_data()
        self.table = self.create_table()
        self.add_save_button()
    
    def create_table_data(self):
        dtype = {
            'Channel.Evidence.Ms1': 'float64', 'Channel.Evidence.Ms2': 'float64',
            'Channel.L': 'float64', 'Channel.M': 'float64',
            'Channel.Q.Value': 'float64', 'Mass.Evidence': 'float64',
            'Ms1.Area': 'float64', 'Ms1.Profile.Corr': 'float64',
            'Ms1.Translated': 'float64', 'Precursor.Normalised': 'float64',
            'Precursor.Quantity': 'float64', 'Precursor.Translated': 'float64',
            'Quantity.Quality': 'float64'
        }
        df = dd.read_csv(os.path.join(self.path, 'report.tsv'), sep='\t', dtype=dtype)
        unique_runs = df['Run'].drop_duplicates().compute()
        return pd.DataFrame({'Run': unique_runs, 'Sample': ['' for _ in unique_runs], 'Treatment': ['' for _ in unique_runs]})
    
    def create_table(self):
        pt = Table(self, dataframe=self.df, showtoolbar=True, showstatusbar=True)
        pt.show()
        return pt
    
    def add_save_button(self):
        save_button = tk.Button(self.main, text='Save', command=self.save_data)
        save_button.pack()
    
    def save_data(self):
        df_to_save = self.table.model.df
        file_path = filedialog.asksaveasfilename(defaultextension='.csv', 
                                                 filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if file_path:
            df_to_save.to_csv(file_path, index=False)
            print(f'Data saved to {file_path}')
    
    def on_closing(self):
        self.main.destroy()