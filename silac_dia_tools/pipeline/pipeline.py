# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 17:53:34 2024

@author: robbi
"""

import os
import pandas as pd
import tkinter as tk
from tkinter import messagebox, filedialog
from pandastable import Table
import json
import dask.dataframe as dd
from icecream import ic

from .utils import manage_directories
from .report import filtering_report, protein_overview_report

from silac_dia_tools.pipeline.preprocessor import Preprocessor 
from silac_dia_tools.pipeline.generate_protein_groups_fix import DiaSis, DynamicDiaSis, DynamicSilac


class Pipeline:
    def __init__(self, path, parameter_file, contains_reference=True, method='dia_sis', pulse_channel="M", meta=None):
        # Assign constructor variables
        self.path = path
        self.parameter_file = parameter_file
        self.pulse_channel = pulse_channel
        self.meta = meta
        self.contains_reference = contains_reference
        self.method = method
   
        # Initialize class variables
        self.relable_with_meta = self._confirm_metadata()
        self.meta_data = self._load_meta_data() if self.relable_with_meta else None
        self.params = self._load_params()
        self.filter_cols = list(self.params['filter_cols'].keys())
       
        # Placeholder variables 
        self.filtered_report = None
        self.contaminants = None
        self.filtered_out_df = None
        
    def _load_params(self):
        json_path = os.path.join(os.path.dirname(__file__), '..', 'configs', self.parameter_file)
        with open(json_path, 'r') as file:
            return json.load(file)

    def _confirm_metadata(self):
        if self.meta is None:
            print("No metadata added, filtering will continue without relabeling")
            return False
        if not isinstance(self.meta, str):
            print("File name is not a string, filtering will continue without relabeling")
            return False
        print("Metadata added, looking for the following file:", self.meta)
        return self._check_directory()

    def _check_directory(self):
        file_list = os.listdir(self.path)
        if self.meta in file_list:
            print(f"CSV file '{self.meta}' found in {self.path}")
            return True
        print(f"CSV file '{self.meta}' not found in the directory.")
        return False

    def _load_meta_data(self):
        return pd.read_csv(os.path.join(self.path, self.meta), sep=',')       
    
    def _save_preprocessing(self):
        manage_directories.create_directory(self.path, 'preprocessing')
        self.contaminants.to_csv(os.path.join(self.path, 'preprocessing', 'contaminants.tsv'), sep='\t')
        self.filtered_out_df.to_csv(os.path.join(self.path, 'preprocessing', 'filtered_out.tsv'), sep='\t')
        self.filtered_report.to_csv(os.path.join(self.path, 'preprocessing', 'filtered_report.tsv'), sep='\t')
        
        # self.LH_df.to_csv(os.path.join(self.path, 'preprocessing', 'light_precursors.tsv'), sep='\t')
        # self.MH_df.to_csv(os.path.join(self.path, 'preprocessing', 'medium_precusors.tsv'), sep='\t')
        # self.href_df.to_csv(os.path.join(self.path, 'preprocessing', 'href.tsv'), sep='\t')
        
    
    def _generate_reports(self):
        # Generate reports for filtering, precursors, and protein groups
        manage_directories.create_directory(self.path, 'reports')
        filtering_report.create_report(self.filtered_report, self.contaminants, self.filtered_out_df, self.path, self.params)
        protein_overview_report.create_report(self.path)
        # precursor_report.create_report(self.formatted_precursors, self.path, self.params, self.method, self.pulse_channel)
        # protein_group_report.create_report(self.protein_groups, self.path, self.params)
        # protein_groups_report_r.create_report(self.path, self.params, self.method)
        print('passsing generate reports steps')
        
    def execute_pipeline(self, generate_report=True):
        self.preprocessor = Preprocessor(self.path, self.params, self.filter_cols, self.contains_reference, self.pulse_channel, self.method, self.meta_data)
        self.filtered_report, self.filtered_out_df, self.contaminants = self.preprocessor.import_report()
        
        if self.method == 'dia_sis':
            self.precursor_rollup = DiaSis(self.path, self.filtered_report)
        elif self.method == 'dynamic_dia_sis':
            self.precursor_rollup = DynamicDiaSis(self.path, self.filtered_report)
        elif self.method == 'dynamic_silac':
            self.precursor_rollup = DynamicSilac(self.path, self.filtered_report, self.pulse_channel)
            # self.generate_report = False
            
        else:
            print('incorrect method')            
        self.precursor_rollup.generate_protein_groups()
        
        self._save_preprocessing()
        
        if generate_report:
            self._generate_reports() 

    def make_metadata(self):
        print("Searching report.tsv for unique runs for metadata, use pop up to enter metadata or copy and past selected runs to a spreadsheet and save as .csv file")
        root = tk.Tk()
        app = TestApp(self.path, root)
        root.protocol("WM_DELETE_WINDOW", app.on_closing)
        app.pack(fill="both", expand=True)  # Ensure the app fills the root window
        root.mainloop()
            
        
class TestApp(tk.Frame):
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

def run_test_app():
    # Example usage
    root = tk.Tk()
    app = TestApp(path='your_path_here', parent=root)
    root.protocol("WM_DELETE_WINDOW", app.on_closing)
    app.pack(fill="both", expand=True)
    root.mainloop()

# if name == "main":
#     run_pipeline() # or run_test_app()

