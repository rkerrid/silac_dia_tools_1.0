from icecream import ic
from silac_dia_tools.workflow.pipeline import Pipeline as pipeline
import pandas as pd 


if __name__ == "__main__": 
   
    
    path = r'G:\My Drive\Data\analysis\20240824 BM pilot 2 vs 3 channel analysis\astral\data\no spike\\'
    # path = r'G:\My Drive\Data\analysis\20240824 BM pilot 2 vs 3 channel analysis\astral\data\spike\\'
    # path = r'C:\phd projects\silac_dia_tools_1.0\tests\unit tests\unit test data\spike\\'
    path = r'G:\My Drive\Data\analysis\20240906 astral ineuron pulse\\'
    
    pipeline = pipeline( f'{path}', 'test_params.json',  method = 'dynamic_silac_dia', pulse_channel="H", metadata_file='meta.csv')

    # pipeline.make_metadata()
    result = pipeline.execute_pipeline()
 
    
 
    