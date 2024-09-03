# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 10:43:55 2024

@author: robbi
"""



class DynamicDiaSis:
    def __init__(self, path, filtered_report):
        self.path = path
        self.filtered_report = filtered_report
        self.update = True
        
        self.formatted_precursors = None
        self.protein_groups = None
        
    def generate_protein_groups(self):
        start_time = time.time()
        
        # formatting SILAC channels
        h_precursors_df, LH_df, MH_df  = self.format_silac_channels(self.filtered_report)
      
        # Calculate global heavy reference df 
        href_df = self.calculate_href_intensities(h_precursors_df)
        
        # calculate protein level ratios
        LH_protein_df = self.compute_protein_level_ratios(LH_df, 'L')
        MH_protein_df = self.compute_protein_level_ratios(MH_df, 'M')
        
        # Merge href with protein groups to generate normalized intensities         
        LH_protein_df, MH_protein_df = self.href_normalization( LH_protein_df, MH_protein_df, href_df)
        
        # oputput data
        self.output_protein_groups( LH_protein_df, MH_protein_df, href_df, self.path)
        
        # save formatted precursors to csv file
        LH_df.to_csv(f'{self.path}formatted_precursors_L.csv', sep=',')
        MH_df.to_csv(f'{self.path}formatted_precursors_M.csv', sep=',')

        end_time = time.time()
        print(f"Time taken to generate protein groups: {end_time - start_time} seconds")
        
    
    def format_silac_channels(self, df):
        # pivot table
        df = df.pivot_table(index=['Run','Protein.Group', 'Precursor.Id'], columns='Label', values = ['Ms1.Translated', 'Precursor.Quantity', 'Precursor.Translated'])
        
        # format href, LH and MH precursors
        href_df = self.format_reference_silac_channels(df)
        LH_df = self.format_LH_silac_channels(df)
        MH_df = self.format_MH_silac_channels(df)
        
        # calculate ratios for LH and MH dfs
        LH_df = self.calculate_precursor_ratios(LH_df, 'L')
        MH_df = self.calculate_precursor_ratios(MH_df, 'M')

        return href_df, LH_df, MH_df    
    
    def format_reference_silac_channels(self, df):
        # Access the 'Ms1.Translated' and 'Precursor.Translated' columns under the 'H' label
        ms1_translated_H = df.loc[:, ('Ms1.Translated', 'H')]
        precursor_translated_H = df.loc[:, ('Precursor.Translated', 'H')]
        
        # Combine into a new DataFrame
        combined_df = pd.DataFrame({
            'Ms1.Translated_H': ms1_translated_H,
            'Precursor.Translated_H': precursor_translated_H
        })
        
        # Reset index to include 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        combined_H_df = combined_df.reset_index()
        
        # Rename columns for clarity
        combined_H_df.columns = ['Run', 'Protein.Group', 'Precursor.Id', 'Ms1.Translated_H', 'Precursor.Translated_H']
        
        # drop rows where there are no valid pairs of Precursor and Ms1 translated
        combined_H_df = combined_H_df.dropna(subset=['Precursor.Translated_H','Ms1.Translated_H'])
        
        return combined_H_df
        
    def format_LH_silac_channels(self, df):
        # Access the 'Ms1.Translated' and 'Precursor.Translated' columns under the 'H' and 'L' label
        ms1_translated_L = df.loc[:, ('Ms1.Translated', 'L')]
        precursor_translated_L = df.loc[:, ('Precursor.Translated', 'L')]
        ms1_translated_H = df.loc[:, ('Ms1.Translated', 'H')]
        precursor_translated_H = df.loc[:, ('Precursor.Translated', 'H')]
        
        # Combine into a new DataFrame
        combined_df = pd.DataFrame({
            'Ms1.Translated_H': ms1_translated_H,
            'Precursor.Translated_H': precursor_translated_H,
            'Ms1.Translated_L': ms1_translated_L,
            'Precursor.Translated_L': precursor_translated_L
        })
        
        # Reset index to include 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        combined_LH_df = combined_df.reset_index()
        
        # Rename columns for clarity
        combined_LH_df.columns = ['Run', 'Protein.Group', 'Precursor.Id', 'Ms1.Translated_H', 'Precursor.Translated_H', 'Ms1.Translated_L', 'Precursor.Translated_L']
               
        return combined_LH_df
    
    def format_MH_silac_channels(self, df):
        # Access the 'Ms1.Translated' and 'Precursor.Translated' columns under the 'H' and 'M' label
        ms1_translated_M = df.loc[:, ('Ms1.Translated', 'M')]
        precursor_translated_M = df.loc[:, ('Precursor.Translated', 'M')]
        ms1_translated_H = df.loc[:, ('Ms1.Translated', 'H')]
        precursor_translated_H = df.loc[:, ('Precursor.Translated', 'H')]
        
        # Combine into a new DataFrame
        combined_df = pd.DataFrame({
            'Ms1.Translated_H': ms1_translated_H,
            'Precursor.Translated_H': precursor_translated_H,
            'Ms1.Translated_M': ms1_translated_M,
            'Precursor.Translated_M': precursor_translated_M
        })
        
        # Reset index to include 'Run', 'Protein.Group', and 'Precursor.Id' as columns
        combined_MH_df = combined_df.reset_index()
        
        # Rename columns for clarity
        combined_MH_df.columns = ['Run', 'Protein.Group', 'Precursor.Id', 'Ms1.Translated_H', 'Precursor.Translated_H', 'Ms1.Translated_M', 'Precursor.Translated_M']
  
        return combined_MH_df
        
    def calculate_precursor_ratios(self, df, channel):
        print('Calculating SILAC ratios based on Ms1.Translated and Precursor.Translated')
        # Calculate ratios for channel        
        df[f'Precursor.Translated {channel}/H'] = df[f'Precursor.Translated_{channel}'] / df['Precursor.Translated_H']
        df[f'Ms1.Translated {channel}/H'] = df[f'Ms1.Translated_{channel}'] / df['Ms1.Translated_H']
        
        columns_to_replace = [f'Precursor.Translated {channel}/H', f'Ms1.Translated {channel}/H']
        df[columns_to_replace] = df[columns_to_replace].replace([np.inf, -np.inf, 0.0], np.nan)
        
        df = df.dropna(subset=[f'Precursor.Translated {channel}/H',f'Ms1.Translated {channel}/H'])
        
        return df
  
    def calculate_href_intensities(self, df):
        
        def combined_median(ms1_series, precursor_series):
            combined_series = np.concatenate([ms1_series, precursor_series])
            combined_series = np.log10(combined_series)  # Log-transform the combined series
            return np.median(combined_series)  # Return the median of the log-transformed values
         
        # Group by protein group and apply the custom aggregation
        grouped = df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
            'href': combined_median(x['Ms1.Translated_H'], x['Precursor.Translated_H']) 
        })).reset_index()
        
        
        return  grouped[['Protein.Group', 'href']]
    
    def compute_protein_level_ratios(self, df, channel):
        runs = df['Run'].unique()
        runs_list = []
        
        for run in tqdm(runs, desc=f'Computing protein level ratios for each run, channel: {channel}'):
            run_df = df[df['Run'] == run]
        
            def combined_median(ms1_series, precursor_series):
                combined_series = np.concatenate([ms1_series, precursor_series])
                # print(combined_series)
                combined_series = np.log10(combined_series)  # Log-transform the combined series
                return np.median(combined_series)  # Return the median of the log-transformed values
             
            # Group by protein group and apply the custom aggregation
            grouped = run_df.groupby(['Protein.Group']).apply(lambda x: pd.Series({
                f'{channel}/H ratio': combined_median(x[f'Ms1.Translated {channel}/H'], x[f'Precursor.Translated {channel}/H'])})).reset_index()
            
            grouped['Run'] = run
            runs_list.append(grouped)
        
        result = pd.concat(runs_list, ignore_index=True)
        
        cols = ['Run','Protein.Group', f'{channel}/H ratio']
         
        return result[cols] 
    
    def href_normalization(self, LH_protein_df, MH_protein_df, href_df):
        # Merge the href_df onto protein groups containing optimized ratios
        merged_df_LH = LH_protein_df.merge(href_df, on='Protein.Group', how='left')
        merged_df_MH = MH_protein_df.merge(href_df, on='Protein.Group', how='left')
        
        # Obtain normalized light intensities by adding the L/H ratio to the heavy refference in log space
        merged_df_LH['L_norm'] = merged_df_LH['L/H ratio'] + merged_df_LH['href']
        merged_df_MH['M_norm'] = merged_df_MH['M/H ratio'] + merged_df_MH['href']
        
        # reverse log data to output protein intensities
        merged_df_LH['L_norm'] = 10**merged_df_LH['L_norm'] 
        merged_df_MH['M_norm'] = 10**merged_df_MH['M_norm']
        href_df['href'] = 10**href_df['href']
        
        return merged_df_LH, merged_df_MH    
    
    def output_protein_groups(self, LH_protein_df, MH_protein_df, href_df, path):
        manage_directories.create_directory(self.path, 'protein_groups')
        print(f'Outputing normalized protein intensities to {path}/protein_groups')
        LH_protein_df = LH_protein_df.rename(columns={'L_norm': 'L'})
        MH_protein_df = MH_protein_df.rename(columns={'M_norm': 'M'})
        
        # Pivoting for 'L'
        l_pivot_df = LH_protein_df.pivot(index='Protein.Group', columns='Run', values='L')
        
        # Pivoting for 'M'
        m_pivot_df = MH_protein_df.pivot(index='Protein.Group', columns='Run', values='M')

        # then output each table to csv for h.href, l.href, m.href
        # h_pivot_df.to_csv(f'{path}/protein_groups/href.csv', sep=',')
        href_df.to_csv(f'{path}/protein_groups/href.csv', sep=',')
        m_pivot_df.to_csv(f'{path}/protein_groups/nsp.csv', sep=',')
        l_pivot_df.to_csv(f'{path}/protein_groups/light.csv', sep=',')