# -*- coding: utf-8 -*-
"""
Created on Sun Jan 14 17:53:34 2024

@author: robbi
"""

import os
import pandas as pd
import json
from icecream import ic
import tkinter as tk

from .utils import manage_directories
from .report import precursor_report

from silac_dia_tools.workflow.preprocessor import Preprocessor 
from silac_dia_tools.workflow.dynamic_dia_sis import DynamicDiaSis
from silac_dia_tools.workflow.dia_sis import DiaSis  
from silac_dia_tools.workflow.dynamic_silac_dia import DynamicSilac
from silac_dia_tools.workflow.meta_data_entry import MetaDataEntry


class Pipeline:
    def __init__(self, path, parameter_file, method='dia_sis', pulse_channel="M", metadata_file=None): # check method input is valid otherwise print method options
        # Assign constructor variables
        self.path = path
        self.parameter_file = parameter_file
        self.pulse_channel = pulse_channel
        self.metadata_file = metadata_file
        self.method = method
   
        # Initialize class variables
        self.relable_with_meta = self._confirm_metadata()
        self.meta_data = self._load_meta_data() if self.relable_with_meta else None
        self.params = self._load_params()
        self.filter_cols = list(self.params['filter_cols'].keys())
       
        # Placeholder variables 
        self.filtered_report = None
        self.contaminants = None
        self.protein_groups = None
        
    def _load_params(self):
        json_path = os.path.join(os.path.dirname(__file__), '..', 'configs', self.parameter_file)
        with open(json_path, 'r') as file:
            return json.load(file)

    def _confirm_metadata(self):
        if self.metadata_file is None:
            print("No metadata added, filtering will continue without relabeling")
            return False
        if not isinstance(self.metadata_file, str):
            print("File name is not a string, filtering will continue without relabeling")
            return False
        print("Metadata added, looking for the following file:", self.metadata_file)
        return self._check_directory()

    def _check_directory(self):
        file_list = os.listdir(self.path)
        if self.metadata_file in file_list:
            print(f"CSV file '{self.metadata_file}' found in {self.path}")
            return True
        print(f"CSV file '{self.metadata_file}' not found in the directory.")
        return False

    def _load_meta_data(self):
        return pd.read_csv(os.path.join(self.path, self.metadata_file), sep=',')       
    
    def _save_preprocessing(self):
        manage_directories.create_directory(self.path, 'preprocessing')
        self.contaminants.to_csv(os.path.join(self.path, 'preprocessing', 'contaminants.csv'), sep=',')
        self.filtered_report.to_csv(os.path.join(self.path, 'preprocessing', 'precursors.csv'), sep=',')
        manage_directories.create_directory(self.path, 'protein_groups')
        self.format_protein_groups(self.protein_groups)     
        self.protein_groups.to_csv(os.path.join(self.path, 'protein_groups', 'protein_groups.csv'), sep=',')
        
    def format_protein_groups(self, protein_groups):
        ''' this function should format protein groups into csv files for light, pulse, heavy?, norm, unnorm, and ratios and save them '''
        return protein_groups
        
    def _generate_reports(self):
        '''this function should call functions from the report module to plot data ralated to IDs, ratios, correlation etc. as well as log data about the run and version etc.'''        
        manage_directories.create_directory(self.path, 'reports')
        precursor_report.create_precursor_report(self.path)
        
    def execute_pipeline(self, generate_report=True):
        self.preprocessor = Preprocessor(self.path,  self.method, self.pulse_channel, self.meta_data)
        self.filtered_report, self.contaminants = self.preprocessor.preprocess()

        if self.method == 'dia_sis':
            self.precursor_rollup = DiaSis(self.path, self.filtered_report)
        elif self.method == 'dynamic_dia_sis':
            self.precursor_rollup = DynamicDiaSis(self.path, self.filtered_report)
        elif self.method == 'dynamic_silac_dia':
            self.precursor_rollup = DynamicSilac(self.path, self.filtered_report)
              
        self.protein_groups = self.precursor_rollup.generate_protein_groups()
     
        self._save_preprocessing()
        self._generate_reports()
        return self.protein_groups
     

    def make_metadata(self):
        print("Searching report.tsv for unique runs for metadata, use pop up to enter metadata or copy and past selected runs to a spreadsheet and save as .csv file")
        root = tk.Tk()
        app = MetaDataEntry(self.path, root)
        root.protocol("WM_DELETE_WINDOW", app.on_closing)
        app.pack(fill="both", expand=True)  # Ensure the app fills the root window
        root.mainloop()


