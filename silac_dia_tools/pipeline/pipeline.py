import os
import dask.dataframe as dd # add to setup
import tkinter as tk # add to setup
from tkinter import messagebox
from tkinter import filedialog
import pandas as pd
from pandastable import Table # add to setup
from icecream import ic 
import json

from silac_dia_tools.pipeline.preprocessor import Preprocessor
from silac_dia_tools.pipeline.silac_formatter import SilacFormatter
from silac_dia_tools.pipeline.precursor_rollup import PrecursorRollup
from silac_dia_tools.pipeline.calculate_intensities import IntensityCalculator
from silac_dia_tools.pipeline.report import filtering_report, precursor_report, protein_group_report, protein_intensities_report
from silac_dia_tools.pipeline.utils import manage_directories


class Pipeline:
    def __init__(self, path, parameter_file, contains_reference = True, pulse_channel="M", meta=None):
        # pipeline variables
        print("Init modified pipeline")
        self.path = path
        self.pulse_channel = pulse_channel
        self.meta = meta
        self.contains_reference = contains_reference
        self.relable_with_meta = self._confirm_metadata()
        self.meta_data = None
        
        if self.relable_with_meta:
            self.meta_data = pd.read_csv(f'{self.path}{self.meta}', sep=',')
            print('Will relabel runs with metadata sample column')
            
        # import parameters
        self.parameter_file = parameter_file
        self.config_dir = os.path.join(os.path.dirname(__file__), '..', 'configs')
        self.params = self._load_params()
        self.filter_cols = list(self.params['apply_filters'].keys()) + ['Lib.PG.Q.Value']
        
        # pipeline objects
        self.preprocessor = Preprocessor(self.path, self.params, self.meta_data)
        self.formatter = SilacFormatter(self.path, self.filter_cols)
        self.precursor_rollup = PrecursorRollup(self.path)
        self.intensity_calculator = IntensityCalculator(self.path, self.contains_reference, self.pulse_channel)
        
        # pipeline outputs
        self.report = None
        self.filtered_report = None
        self.filtered_out = None
        self.contaminants = None
        
        self.formatted_precursors = None
        self.filtered_precursors = None
        self.protein_groups = None
        
        self.unnormalized_total_intensities = None
        self.unnormalized_nsp_intensities = None
        self.light_unormalized_intensities = None
        
        self.href_total_intensities = None
        self.href_nsp_intensities = None
        self.href_nsp_light = None
        
        self.dlfq_nsp_intensities = None
        self.dlfq_total_intensities = None
        self.dlfq_total_light = None
    
    def _load_params(self):
        json_path = os.path.join(self.config_dir, self.parameter_file)
        with open(json_path, 'r') as f:
            params = json.load(f)
            return params 
    
    def _confirm_metadata(self):
        if self.meta is None:
            print("No metadata added, filtering will continue without relabeling")
            return False
        elif isinstance(self.meta, str):
            print("Metadata added, looking for the following file:", self.meta)
            meta_exists = self._check_directory()
            if meta_exists:
                return True
            else:
                print(f"Cannot find {self.meta} in {self.path}, check file name and location") 
                print("Filtering will continue without relabeling")
        else:
            print("File name is not a string, filering will continue without relabeling")
            return False
    
    def _check_directory(self):
        file_list = os.listdir(self.path)
        # Iterate through the list of filenames and check for a match
        found = False
        for filename in file_list:
            if filename == self.meta:
                found = True
                print(f"CSV file '{self.meta}' found in {self.path}")
                return True
        if not found:
            print(f"CSV file '{self.meta}' not found in the directory.")
            return False
        
    def preprocess_href(self):
        self.report = self.preprocessor.import_report()
        print('finished reading in file')
        self.formatted_precursors = self.formatter.format_silac_channels(self.report)
        self.filtered_report, self.contaminants, self.filtered_out = self.preprocessor.filter_formatted(self.formatted_precursors)
        
        self.filtered_precursors  = self.filtered_report[self.filtered_report['quantity type'] == 'Ms1.Translated']
        self.contaminants  = self.contaminants[self.contaminants['quantity type'] == 'Ms1.Translated']
        self.filtered_out = self.filtered_out[self.filtered_out['quantity type'] == 'Ms1.Translated']
        
        self._roll_up_to_protein_level()
        self.save_preprocessing()
        self._output_unnormalized()
       
        self._output_href()
        
    def preprocess_dlfq(self):
        self.report = self.preprocessor.import_report()
        # self.report = pd.read_csv(f'G:/My Drive/Data/data/testing pipeline dev/eif4f/subset_report.csv', sep='\t')
        print('finished reading in file')
        self.formatted_precursors = self.formatter.format_silac_channels(self.report)
        self.filtered_report, self.contaminants, self.filtered_out = self.preprocessor.filter_formatted(self.formatted_precursors)
        
        self.filtered_precursors  = self.filtered_report[self.filtered_report['quantity type'] == 'Ms1.Translated']
        self.contaminants  = self.contaminants[self.contaminants['quantity type'] == 'Ms1.Translated']
        self.filtered_out = self.filtered_out[self.filtered_out['quantity type'] == 'Ms1.Translated']
        
        self._roll_up_to_protein_level()
        self.save_preprocessing()
        self._output_unnormalized()
       
        self._output_dlfq()
    
    def make_metadata(self):
        root = tk.Tk()
        app = TestApp(self.path, root)
        root.protocol("WM_DELETE_WINDOW", app.on_closing)
        app.pack(fill="both", expand=True)  # Ensure the app fills the root window
        root.mainloop()
        
    def ask_user(self):
        print('ask user')
        root = tk.Tk()
        root.withdraw()
        
        response = messagebox.askyesno('Choose Option', 'Do you want to create a metadata file?')
        
        if response:
            print('Add metadata and save file as meta.csv in report.tsv directory')
            self.meta = f'{self.path}meta.csv'
            self.make_metadata()
         
        else:
            print('Continue without creating metadata file')
            
    def _preprocess(self):
        self.filtered_report, self.contaminants, self.filtered_out = self.preprocessor.import_and_filter()
        
    def _format_channels(self):
        self.formatted_precursors = self.formatter.format_silac_channels(self.filtered_report)
    
    def _roll_up_to_protein_level(self):
        print('uses filtered precursors rather than report')
        self.protein_groups = self.precursor_rollup.calculate_protein_level_ratios(self.filtered_precursors) #jsut changed this to filtered rreport rather than precursors dut to bug?
        
        
    def _output_unnormalized(self):   
       self.unnormalized_total_intensities, self.unnormalized_nsp_intensities, self.light_unormalized_intensities = self.intensity_calculator.output_unnorm(self.protein_groups)
       
    def _output_href(self):
       self.href_total_intensities, self.href_nsp_intensities, self.href_nsp_light = self.intensity_calculator.output_href(self.protein_groups)
       
    def _output_dlfq(self):
       self.dlfq_total_intensities, self.dlfq_nsp_intensities, self.dlfq_total_light = self.intensity_calculator.output_dlfq(self.protein_groups, self.filtered_precursors)

    def save_preprocessing(self):
        manage_directories.create_directory(self.path, 'preprocessing')
        
        print('Saving silac_precursors.tsv')
        self.filtered_precursors.to_csv(f'{self.path}preprocessing/silac_precursors.tsv',sep='\t')
        print('Saving protein_ratios.csv')
        self.filtered_report.to_csv(f'{self.path}preprocessing/protein_ratios.csv',sep='\t')
       
        
    def run_href_pipeline(self):
        self._preprocess()
        self._format_channels()
        self._roll_up_to_protein_level()
        self.save_preprocessing()
        self._output_unnormalized()
       
        self._output_href()
        
    def run_dlfq_pipeline(self):
        self._preprocess()
        self._format_channels()
        self._roll_up_to_protein_level()
        self.save_preprocessing()
        self._output_unnormalized()
        
        self._output_dlfq()
        
    def generate_reports(self): 
        print('Beginning filtering report')
        filtering_report.create_report(self.filtered_report, self.contaminants, self.filtered_out, self.path, self.params)
        print('Beginning precursor report')
        precursor_report.create_report(self.formatted_precursors, self.path, self.params)
        print('Beginning protein group report')
        protein_group_report.create_report(self.protein_groups, self.path, self.params)
        # print('Beginning protein intensities report')
        # file_list = [f for f in os.listdir(f'{self.path}protein intensities') if os.path.isfile(os.path.join(f'{self.path}protein intensities', f))]
        if not self.contains_reference:
            protein_intensities_report.create_report(self.path, self.params)
        

class TestApp(tk.Frame):
    def __init__(self, path, parent=None):
        self.path = path
        self.parent = parent
        tk.Frame.__init__(self, parent)
        self.main = self.master
        self.main.geometry('600x400+200+100')
        self.main.title('Data Entry Table')
        self.df = self.create_table_data(path)
        self.table = self.create_table()

    def create_table_data(self, path):
        # Predefined 'run' list
        print('Beginning dd import')
        dtype={'Channel.Evidence.Ms1': 'float64',
                'Channel.Evidence.Ms2': 'float64',
                'Channel.L': 'float64',
                'Channel.M': 'float64',
                'Channel.Q.Value': 'float64',
                'Mass.Evidence': 'float64',
                'Ms1.Area': 'float64',
                'Ms1.Profile.Corr': 'float64',
                'Ms1.Translated': 'float64',
                'Precursor.Normalised': 'float64',
                'Precursor.Quantity': 'float64',
                'Precursor.Translated': 'float64',
                'Quantity.Quality': 'float64'}
        
        df = dd.read_csv(f'{path}report.tsv', sep='\t',dtype=dtype)
        unique_runs = df['Run'].drop_duplicates().compute()
        print(f'unique runs: {unique_runs}')
        # Create DataFrame
        df = pd.DataFrame({'Run': unique_runs, 'Sample': ['' for _ in unique_runs], 'Treatment': ['' for _ in unique_runs]})
        return df

    def create_table(self):
        pt = Table(self, dataframe=self.df, showtoolbar=True, showstatusbar=True)
        pt.show()
        return pt
    
    def add_save_button(self):
        save_button = tk.Button(self.main, text='Save', command=self.save_data)
        save_button.pack()

    def save_data(self):
        # Get DataFrame from table
        df_to_save = self.table.model.df

        # Open file dialog to select path and file name
        file_path = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if file_path:  # Check if a file path was selected
            df_to_save.to_csv(file_path, index=False)
            print(f'Data saved to {file_path}')

    def on_closing(self):
        # This function is called when the window is closed
        # You can add code here to handle the DataFrame
        print(self.table.model.df)  # Example: Print the DataFrame
        self.main.destroy()




