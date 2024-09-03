from icecream import ic
from silac_dia_tools.workflow.pipeline import Pipeline as pipeline
import pandas as pd 


if __name__ == "__main__": 
   
    
    path = 'G:/My Drive/Data/data/20240903 test pipeline with tripple bm/n/'
    
    pipeline = pipeline( f'{path}', 'test_params.json',  method = 'dynamic_silac_dia', pulse_channel="M", metadata_file='meta.csv')

    # pipeline.make_metadata()
    result = pipeline.execute_pipeline()
 
    
 
    