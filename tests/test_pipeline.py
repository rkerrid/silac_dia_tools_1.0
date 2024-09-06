from icecream import ic
from silac_dia_tools.workflow.pipeline import Pipeline as pipeline
import pandas as pd 


if __name__ == "__main__": 
   
    
    # path = r'G:\My Drive\Data\analysis\20240824 BM pilot 2 vs 3 channel analysis\astral\data\no spike\\'
    path = r'G:\My Drive\Data\analysis\20240824 BM pilot 2 vs 3 channel analysis\astral\data\spike\\'
    # path = 'G:/My Drive/Data/data/20240804 jose neurones 48 silac astral/'
    
    pipeline = pipeline( f'{path}', 'test_params.json',  method = 'dynamic_dia_sis', pulse_channel="M", metadata_file='meta.csv')

    # pipeline.make_metadata()
    result = pipeline.execute_pipeline()
 
    
 
    