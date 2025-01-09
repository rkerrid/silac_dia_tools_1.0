from icecream import ic
from silac_dia_tools.workflow.pipeline import Pipeline as pipeline
import pandas as pd 


if __name__ == "__main__": 
   
     
   path = r'G:\My Drive\Data\data\20250109 bm20 nospike dian192\report.tsv\\'
    
   pipeline = pipeline( f'{path}',  method = 'dynamic_silac_dia', pulse_channel="M", metadata_file='meta.csv')
    
   pipeline.execute_pipeline()