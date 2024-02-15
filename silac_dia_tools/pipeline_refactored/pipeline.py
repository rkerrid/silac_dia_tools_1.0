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

from .utils import manage_directories
from .report import filtering_report, precursor_report, protein_group_report, protein_intensities_report

from silac_dia_tools.pipeline_refactored.preprocessor import Preprocessor 
from silac_dia_tools.pipeline_refactored.generate_protein_groups import HrefRollUp, LfqRollUp 
# to do below
# from silac_dia_tools.pipeline_refactored.silac_formatter import SilacFormatter 
# from silac_dia_tools.pipeline_refactored.calculate_intensities import IntensityCalculator
# from silac_dia_tools.pipeline_refactored.precursor_rollup import PrecursorRollup

# Import custom modules
# from silac_dia_tools.pipeline_refactored import (
#     preprocessor, silac_formatter, precursor_rollup,
#     calculate_intensities
# )


class Pipeline:
    def __init__(self, path, parameter_file, contains_reference=True, pulse_channel="M", meta=None):
        # Assign constructor variables
        self.path = path
        self.parameter_file = parameter_file
        self.pulse_channel = pulse_channel
        self.meta = meta
        self.contains_reference = contains_reference
        
        # Initialize class variables
        self.relable_with_meta = self._confirm_metadata()
        self.meta_data = self._load_meta_data() if self.relable_with_meta else None
        self.params = self._load_params()
        self.filter_cols = list(self.params['filter_cols'].keys())
       
        # Placeholder variables 
        self.precursor_rollup = None
        
        self.filtered_report = None
        self.contaminants = None
        self.filtered_out_df = None
        self.formatted_precursors = None
        self.protein_groups = None
        
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
        self.filtered_report.to_csv(os.path.join(self.path, 'preprocessing', 'filtered_report.tsv'), sep='\t')
        self.formatted_precursors.to_csv(os.path.join(self.path, 'preprocessing', 'silac_precursors.tsv'), sep='\t')
        self.protein_groups.to_csv(os.path.join(self.path, 'preprocessing', 'protein_groups.csv'), sep=',')
    
    def _generate_reports(self):
        # # For the report since we are just using precursor translated these dfs need to be subsetted first
        # filter_type = 'Precursor.Translated'
        # self.filtered_precursors = self.filtered_report[self.filtered_report['quantity type'] == filter_type]
        # self.contaminants = self.contaminants[self.contaminants['quantity type'] == filter_type]
        # self.filtered_out = self.filtered_out[self.filtered_out['quantity type'] == filter_type]
        
        # Then call the reporting modules for each preprocessing step
        filtering_report.create_report(self.filtered_report, self.contaminants, self.filtered_out_df, self.path, self.params)
        precursor_report.create_report(self.formatted_precursors, self.path, self.params)
        protein_group_report.create_report(self.protein_groups, self.path, self.params)

    def execute_pipeline(self, generate_report=True):
        self.preprocessor = Preprocessor(self.path, self.params, self.filter_cols, self.contains_reference, self.pulse_channel, self.meta_data)
        self.filtered_report, self.filtered_out_df, self.contaminants = self.preprocessor.import_report()
        
        if self.contains_reference:
            self.precursor_rollup = HrefRollUp(self.path, self.filtered_report)
        else:
            self.precursor_rollup = LfqRollUp(self.path, self.filtered_report, self.pulse_channel)
            
        self.formatted_precursors, self.protein_groups = self.precursor_rollup.generate_protein_groups()
        
        if generate_report:
            self._save_preprocessing()
            self._generate_reports() 
          
          
class TestApp(tk.Frame):
    def init(self, path, parent=None):
        super().init(parent)
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


def run_pipeline():
    # Example usage
    pipeline = Pipeline(path='your_path_here', parameter_file='your_parameter_file.json')
    pipeline.preprocess_pipeline(method='href')
    pipeline.generate_reports()

def run_test_app():
    # Example usage
    root = tk.Tk()
    app = TestApp(path='your_path_here', parent=root)
    root.protocol("WM_DELETE_WINDOW", app.on_closing)
    app.pack(fill="both", expand=True)
    root.mainloop()

# if name == "main":
#     run_pipeline() # or run_test_app()

