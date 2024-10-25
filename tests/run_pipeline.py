from icecream import ic
from silac_dia_tools.workflow.pipeline import Pipeline as pipeline
import pandas as pd 


if __name__ == "__main__": 
   
    
    # # path = r'G:\My Drive\Data\analysis\20240824 BM pilot 2 vs 3 channel analysis\astral\data\no spike\\'
    # # path = r'G:\My Drive\Data\analysis\20240824 BM pilot 2 vs 3 channel analysis\astral\data\spike\\'
    # path = r'G:\My Drive\Data\main experiments\20240912 Tripple bm figures\data\20240912 exploris new version\spike\\'
    # # path = r'G:\My Drive\Data\main experiments\20240912 Tripple bm figures\data\20240912 astral new version\spike\\'
    
    # pipeline = pipeline( f'{path}',  method = 'dynamic_dia_sis', pulse_channel="M", metadata_file='meta.csv')

    # # pipeline.make_metadata()
    # result = pipeline.execute_pipeline()
 
    
    # path = r'G:\My Drive\Data\main experiments\20241025 astral BM\BM38 nospike\\'
    # path = r'G:\My Drive\Data\main experiments\20241025 astral BM\BM23 nospike\\'
    # path = r'G:\My Drive\Data\main experiments\20241025 astral BM\BM23 spike\\'
    # path = r'G:\My Drive\Data\main experiments\20241017 astral dilution BM\nospike\\'
    path = r'G:\My Drive\Data\main experiments\20241024 HCT116 pulse poc\\'
    
    pipeline = pipeline( f'{path}',  method = 'dynamic_dia_sis', pulse_channel="M", metadata_file='meta.csv')
    # pipeline.make_metadata()
    result = pipeline.execute_pipeline()