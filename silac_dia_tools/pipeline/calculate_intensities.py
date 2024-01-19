import pandas as pd
from silac_dia_tools.pipeline.utils import dlfq_functions as dlfq
from silac_dia_tools.pipeline.utils import manage_directories


class IntensityCalculator:
    def __init__(self, path, contains_reference, pulse_channel):
        self.path = path
        self.pulse_channel = pulse_channel
        self.contains_reference = contains_reference

    def output_href(self, ratios):
        manage_directories.create_directory(self.path,'protein intensities')
        print('Calculating href intensities')
     
        # Generate href df
        h_ref = ratios.groupby('Protein.Group')['H intensity'].median()
        h_ref = h_ref.reset_index().rename(columns={'H intensity': 'h_ref'})
        
        # merge href onto ratios
        merged_df_h = ratios.merge(h_ref, on='Protein.Group', how='inner')

        merged_df_h['H normalized total intensity'] = merged_df_h['h_ref'] / merged_df_h['H to stack ratio']
        merged_df_h['H intensity'] = merged_df_h['H normalized total intensity'] * merged_df_h['H to stack ratio']
        merged_df_h['M intensity'] = merged_df_h['H normalized total intensity'] * merged_df_h['M to stack ratio']
        merged_df_h['L intensity'] = merged_df_h['H normalized total intensity'] * merged_df_h['L to stack ratio']
            
        # Assign intensities to relevant columns
        merged_df_h['Total intensity'] =  merged_df_h['L intensity'] + merged_df_h['M intensity']
        merged_df_h['NSP intensity'] = merged_df_h['M intensity']
        
        # Generate subsetted dfs based on channels
        light_hnorm =  merged_df_h[['Run', 'Protein.Group', 'L intensity']]
        nsp_hnorm = merged_df_h[['Run', 'Protein.Group', 'NSP intensity']]
        total_hnorm = merged_df_h[['Run', 'Protein.Group', 'Total intensity']]
        reference = merged_df_h[['Run', 'Protein.Group', 'H intensity']]
        
        # Pivot tables to output format
        light_hnorm = light_hnorm.pivot(index='Protein.Group', columns='Run', values='L intensity')
        nsp_hnorm = nsp_hnorm.pivot(index='Protein.Group', columns='Run', values='NSP intensity')
        total_hnorm = total_hnorm.pivot(index='Protein.Group', columns='Run', values='Total intensity')
        reference = reference.pivot(index='Protein.Group', columns='Run', values='H intensity')
        
        light_hnorm.to_csv(f'{self.path}protein intensities/light_href.csv', sep=',')  
        nsp_hnorm.to_csv(f'{self.path}protein intensities/nsp_href.csv', sep=',')
        total_hnorm.to_csv(f'{self.path}protein intensities/total_href.csv', sep=',')
        reference.to_csv(f'{self.path}protein intensities/reference_href.csv', sep=',')
    
        print('Saved H reference normalized protein intensities')
        return total_hnorm, nsp_hnorm, light_hnorm

                
    def output_unnorm(self, ratios, contains_reference=True):
        # self.create_protein_intensity_directory()
        manage_directories.create_directory(self.path,'protein intensities')
        print('Calculating unnormalized intensities')
   
        # Assign intensities to relevant columns
        nsp_channel = self.pulse_channel + ' intensity'
        ratios['Total intensity'] = ratios['L intensity'] +  ratios[nsp_channel]
        ratios['NSP intensity'] = ratios[nsp_channel]
        
        # Generate subsetted dfs based on channels
        light_unnorm = ratios[['Run', 'Protein.Group', 'L intensity']]
        nsp_unnorm = ratios[['Run', 'Protein.Group', 'NSP intensity']]
        total_unnorm = ratios[['Run', 'Protein.Group', 'Total intensity']]
        
        # Pivot tables to output format
        light_unnorm = light_unnorm.pivot(index='Protein.Group', columns='Run', values='L intensity')
        nsp_unnorm = nsp_unnorm.pivot(index='Protein.Group', columns='Run', values='NSP intensity')
        total_unnorm = total_unnorm.pivot(index='Protein.Group', columns='Run', values='Total intensity')
        
        # Save tables to CSV
        light_unnorm.to_csv(f'{self.path}protein intensities/light_unnorm.csv', sep=',')
        nsp_unnorm.to_csv(f'{self.path}protein intensities/nsp_unnorm.csv', sep=',')
        total_unnorm.to_csv(f'{self.path}protein intensities/total_unnorm.csv', sep=',')
        if contains_reference:
            reference = ratios[['Run', 'Protein.Group', 'H intensity']]
            reference = reference.pivot(index='Protein.Group', columns='Run', values='H intensity')
            reference.to_csv(f'{self.path}protein intensities/reference_unnorm.csv', sep=',')
        print('Saved unnormalized protein intensities')

        return nsp_unnorm, total_unnorm, light_unnorm

    def output_dlfq(self, ratios, silac_precursors):
        manage_directories.create_directory(self.path,'protein intensities')
        print('Calculating dlfq intensities')
        nsp_ratio = self.pulse_channel + ' to stack ratio'
        
        silac_precursors = silac_precursors[silac_precursors['quantity type']== 'Ms1.Translated']
        silac_precursors.to_csv(f'{self.path}preprocessing/silac_precursors_dlfq_in.tsv', sep='\t')
        dlfq.run_lfq( f'{self.path}preprocessing/silac_precursors_dlfq_in.tsv', 
                  file =  f'{self.path}preprocessing/dlfq_protein_intensities.tsv',
                  num_cores=1)
        
        # After directLFQ finishes running, read in results, format, and merge onto ratios df
        lfq_df = pd.read_csv(f'{self.path}preprocessing/dlfq_protein_intensities.tsv', sep='\t')
        lfq_df = lfq_df.melt(id_vars=['protein'], var_name = 'Run', value_name = 'Intensity')
        lfq_df.rename(columns = {'protein': 'Protein.Group'}, inplace=True)
        merged_df = ratios.merge(lfq_df, on=['Protein.Group','Run'], how='inner')
        
        # Assign intensities to relevant columns
        merged_df['L intensity'] = merged_df['L to stack ratio'] * merged_df['Intensity']
        merged_df['Total intensity'] = (merged_df['L to stack ratio'] * merged_df['Intensity']) + (merged_df[nsp_ratio] * merged_df['Intensity'])
        merged_df['NSP intensity'] = (merged_df[nsp_ratio] * merged_df['Intensity'])
        
        # Generate subsetted dfs based on channels
        total_lfq = merged_df[['Run', 'Protein.Group', 'Total intensity']]
        nsp_lfq = merged_df[['Run', 'Protein.Group', 'NSP intensity']]
        light_lfq = merged_df[['Run', 'Protein.Group', 'L intensity']]
        
        # Pivot tables to output format
        total_lfq = total_lfq.pivot(index='Protein.Group', columns='Run', values='Total intensity')
        nsp_lfq = nsp_lfq.pivot(index='Protein.Group', columns='Run', values='NSP intensity')
        light_lfq = light_lfq.pivot(index='Protein.Group', columns='Run', values='L intensity')
        
        # nsp to light ratio and light to H ratio
        
        # Save tables to CSV
        total_lfq.to_csv(f'{self.path}protein intensities/total_dlfq.csv', sep=',')
        nsp_lfq.to_csv(f'{self.path}protein intensities/nsp_dlfq.csv', sep=',')
        light_lfq.to_csv(f'{self.path}protein intensities/light_dlfq.csv', sep=',')
        print('Saved directLFQ normalized protein intensities')

        return total_lfq, nsp_lfq, light_lfq


