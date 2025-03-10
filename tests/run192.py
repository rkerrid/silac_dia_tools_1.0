from icecream import ic
from silac_dia_tools.workflow_192.pipeline import Pipeline as pipeline
import pandas as pd 


if __name__ == "__main__": 
   
     
   path = r'G:\My Drive\Data\PhD results chapters\Three channel BM on all machines\data\20250131 BM astral 40spdwz diann 92\\'
   path = path = r'W:\CGI\20250206_Popeye_DIA_T2\DIANN_1dot8dot2_output\\'
   # path = r'G:\My Drive\Data\data\20240909 BM pilot astral\no spike\\'
   pipeline = pipeline( f'{path}',  method = 'dynamic_silac_dia', pulse_channel="H", metadata_file='meta.csv')
    
   pipeline.execute_pipeline()
   # pipeline.make_metadata()
   
   
   ['File.Name', 
    'Run', 
    'Protein.Group', 
    'Protein.Ids', 
    'Protein.Names', 
    'Genes', 
    'PG.Quantity', 
    'PG.Normalised', 
    'PG.MaxLFQ', 
    'Genes.Quantity',
    'Genes.Normalised', 
    'Genes.MaxLFQ', 
    'Genes.MaxLFQ.Unique',
    'Modified.Sequence', 
    'Stripped.Sequence', 
    'Precursor.Id', 
    'Precursor.Charge',
    'Q.Value', 
    'PEP', 
    'Global.Q.Value',
    'Protein.Q.Value', 
    'PG.Q.Value', 
    'Global.PG.Q.Value',
    'GG.Q.Value', 
    'Translated.Q.Value', 
    'Proteotypic', 
    'Precursor.Quantity',
    'Precursor.Normalised',
    'Quantity.Quality', 
    'RT', 'RT.Start',
    'RT.Stop', 'iRT', 
    'Predicted.RT', 'Predicted.iRT',
    'First.Protein.Description', 
    'Lib.Q.Value', 'Lib.PG.Q.Value',
    'Ms1.Profile.Corr', 
    'Ms1.Area', 'Ms1.Normalised', 
    'Normalisation.Factor', 
    'Evidence', 'Spectrum.Similarity', 
    'Averagine', 'Mass.Evidence', 'CScore', 
    'Fragment.Quant.Raw', 'Fragment.Correlations', 
    'MS2.Scan', 'Channel.Evidence.Ms1', 'Channel.Evidence.Ms2', 
    'Channel.Q.Value', 'Channel.L', 'Channel.M', 
    'Channel.decoy', 'Precursor.Mz', 'Fragment.Info', 
    'Lib.Index', 'IM', 'iIM', 'Predicted.IM', 'Predicted.iIM']