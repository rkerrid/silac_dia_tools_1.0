from icecream import ic
from silac_dia_tools.pipeline.pipeline import Pipeline as pipeline
import pandas as pd 


if __name__ == "__main__": 
   
    
    path = 'G:/My Drive/Data/data/20240624 starvation pilot/'
    path = 'G:/My Drive/Data/main experiments/20240712 starvation pilot/'
    
    pipeline = pipeline( f'{path}', 'test_params.json', contains_reference = True, method = 'dynamic_dia_sis', pulse_channel="M", meta='meta.csv')

    # pipeline.make_metadata()
    pipeline.execute_pipeline()
 
    
 
    