from icecream import ic
from silac_dia_tools.workflow.pipeline import Pipeline as pipeline
import pandas as pd 


if __name__ == "__main__": 
   
    
    path = 'G:/My Drive/Data/data/20240903 test pipeline with tripple bm/n/'
    path = 'G:/My Drive/Data/data/20240804 jose neurones 48 silac astral/'
    
    pipeline = pipeline( f'{path}', 'test_params.json',  method = 'dynamic_silac_dia', pulse_channel="H", metadata_file='meta.csv')

    # pipeline.make_metadata()
    result = pipeline.execute_pipeline()
 
    
 
    