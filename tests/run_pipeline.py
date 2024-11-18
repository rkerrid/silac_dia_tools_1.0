from icecream import ic
from silac_dia_tools.workflow.pipeline import Pipeline as pipeline
import pandas as pd 


if __name__ == "__main__": 
   
     

    

    # path = r'G:\My Drive\Data\main experiments\20241012 timsTOF BM\SPD30 dil spike\\'
    # pipeline = pipeline( f'{path}',  method = 'dynamic_dia_sis', pulse_channel="M", metadata_file='meta.csv')
    # result = pipeline.execute_pipeline()
    
    # path = r'G:\My Drive\Data\main experiments\20241012 timsTOF BM\SPD30 spike\\'
    # pipeline = pipeline( f'{path}',  method = 'dynamic_dia_sis', pulse_channel="M", metadata_file='meta.csv')
    # result = pipeline.execute_pipeline()
    
    path = r'G:\My Drive\Data\main experiments\20241017 astral dilution BM\nospike\\'
    path = r'W:\RJK\Flo_20241118\\'
    pipeline = pipeline( f'{path}',  method = 'dynamic_silac_dia', pulse_channel="H", metadata_file='meta.csv')
    result = pipeline.execute_pipeline()
     
    path = r'W:\RJK\Flo_20241118\\'

    # path = r'G:\My Drive\Data\main experiments\20241012 timsTOF BM\SPD30 dil nospike\\'
    # pipeline = pipeline( f'{path}',  method = 'dynamic_silac_dia', pulse_channel="M", metadata_file='meta.csv')
    # result = pipeline.execute_pipeline()
    
    # path = r'G:\My Drive\Data\main experiments\20241012 timsTOF BM\SPD30 nospike\\'
    # pipeline = pipeline( f'{path}',  method = 'dynamic_silac_dia', pulse_channel="M", metadata_file='meta.csv')
    # result = pipeline.execute_pipeline()
    
    # path = r'G:\My Drive\Data\main experiments\20241012 timsTOF BM\SPD60 nospike\\'
    # pipeline = pipeline( f'{path}',  method = 'dynamic_silac_dia', pulse_channel="M", metadata_file='meta.csv')
    # result = pipeline.execute_pipeline()