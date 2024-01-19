# # -*- coding: utf-8 -*-
# """
# Created on Wed Nov 29 11:52:35 2023

# @author: rkerrid
# """

from icecream import ic
from silac_dia_tools.pipeline.pipeline import Pipeline as pileline
import pandas as pd 

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
if __name__ == "__main__":
    
 
   
    
    path = 'G:/My Drive/Data/data/eif4g optimization/'
    # path = 'G:/My Drive/Data\data/240112 poc4 test/new pipeline and stats/'
    pipeline = pileline( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta='meta.csv')
    # pipeline.make_metadata()
    pipeline.preprocess_href() # in href mode
    pipeline.generate_reports()
    
    # path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/new h filter test/' 
    # # path = 'G:/My Drive/Data\data/240112 poc4 test/new pipeline and stats/'
    # pipeline = pileline( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta='meta.csv')
    # # pipeline.make_metadata()
    # pipeline.preprocess_href() # in href mode
    # pipeline.generate_reports()
    
    
    # path = 'G:/My Drive/Data/data/240112 poc4 test/new pipeline new stats/N/'
    # pipeline = pileline( f'{path}', 'filtering_parameters_strict.json', contains_reference = False, pulse_channel="M", meta='meta.csv')
    # # pipeline.make_metadata()
    # pipeline.preprocess_dlfq() # in href mode
    # pipeline.generate_reports()
    
    path = 'G:/My Drive/Data/data/eIF4F optimization/'
    # path = 'G:/My Drive/Data\data/240112 poc4 test/new pipeline and stats/'
    pipeline = pileline( f'{path}', 'filtering_parameters_strict.json', contains_reference = True, pulse_channel="M", meta='meta.csv')
    # pipeline.make_metadata()
    pipeline.preprocess_href() # in href mode
    pipeline.generate_reports()
    

# from icecream import ic
# from silac_dia_tools.pipeline.refactored.pipeline_refactored import Pipeline 
# import pandas as pd 

# data_path = 'G:/My Drive/Data/data/test_new_pipeline/H/'
# parameter_file_name = 'filtering_parameters_strict.json'

# def run_pipeline():
#     pipeline = Pipeline(path=data_path, parameter_file=parameter_file_name)
#     pipeline.preprocess_pipeline(method='href')  # You can change 'href' to 'dlfq' based on your need
#     # pipeline.generate_reports()

# # def run_test_app():
# #     root = tk.Tk()
# #     app = TestApp(path=data_path, parent=root)
# #     root.protocol("WM_DELETE_WINDOW", app.on_closing)
# #     app.pack(fill="both", expand=True)
# #     root.mainloop()

# if __name__ == "__main__":
#     run_pipeline()  # Comment this out if you want to run the TestApp GUI instead
#     # run_test_app()

    
   
   
