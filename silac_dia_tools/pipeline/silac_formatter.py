
import pandas as pd
import numpy as np
from icecream import ic


class SilacFormatter:
    def __init__(self, path, filter_cols):
        self.path = path
        self.report = None
        self.parsed_df = None
        self.combined_precursors = None
        self.stacked_intensities = None
        self.filter_cols = filter_cols
        
    def format_silac_channels(self, report):
        print('Beggining formatiing silac channels')
        # print('before parsing')
        self.report = report
        # ic(self.report)
        self.parsed_df = self.parse_data_for_channel_info(self.report, self.filter_cols)
       
        self.combined_precursors = self.combine_modified_precursors()
        self.stacked_intensities = self.stack_intensities()
        
        print('Finished formatiing silac channels')
        return self.stacked_intensities

    def parse_data_for_channel_info(self, report, filter_cols):
        print('Parsing data for SILAC intensities')
        report['Label'] = report['Precursor.Id'].str.extract(r'\(SILAC-(K|R)-([HML])\)')[1]    
        report['Precursor'] = report['Stripped.Sequence'].astype(str) + report['Precursor.Charge'].astype(str)
        # Splitting the data into two, one for 'Ms1.Translated' and another for 'Precursor.Translated'
        ms1_translated_df = report.copy()
        ms1_translated_df['intensity'] = report['Ms1.Translated']
        ms1_translated_df['quantity type'] = 'Ms1.Translated'
    
        precursor_translated_df = report.copy()
        precursor_translated_df['intensity'] = report['Precursor.Translated']
        precursor_translated_df['quantity type'] = 'Precursor.Translated'
        
        # Concatenate the two dataframes to create the parsed dataframe
        parsed_df = pd.concat([ms1_translated_df, precursor_translated_df], ignore_index=True)
        return parsed_df[['Run', 'Protein.Group', 'Precursor.Id', 'Precursor', 'Label', 'intensity', 'quantity type', 'Precursor.Quantity'] + filter_cols]

    # def combine_modified_precursors(self):
    #     print('combining modified precursors')
    #     # Define aggregation functions for columns
    #     keys = self.filter_cols
    #     agg_functions = {key: 'first' for key in keys}
    #     agg_functions['intensity'] = 'first'
    #     # Aggregate data using groupby and agg function
    #     agg_df = self.parsed_df.groupby(['Run', 'Protein.Group',  'Precursor', 'quantity type', 'Label']).agg(agg_functions).reset_index()

    #     # Pivot the table to get 'H', 'M', 'L' intensities in separate columns
    #     pivoted_df = agg_df.pivot_table(index=['Run', 'Protein.Group', 'Precursor', 'quantity type'], columns='Label', values='intensity', fill_value=0).reset_index()

    #     # Rename columns based on the label
    #     pivoted_df.columns.name = None
    #     pivoted_df = pivoted_df.rename(columns={'H': 'H intensity', 'M': 'M intensity', 'L': 'L intensity'})

    #     # Merge the pivoted_df with the original df to get other columns back
    #     df = pd.merge(agg_df.drop(columns=['Label', 'intensity']), pivoted_df, on=['Run', 'Protein.Group',  'Precursor', 'quantity type'])
        
    #     # Assuming df is your DataFrame
    #     # 'column_to_keep_highest' is the column for which you want to keep the highest value
    #     # 'subset_columns' is a list of columns based on which you want to drop duplicates
        
    #     df = df.sort_values(self.filter_cols, ascending=False).drop_duplicates(subset=['Run', 'Protein.Group',  'Precursor', 'quantity type'])

    #     df = df.drop_duplicates()
    #     return df
    
    def combine_modified_precursors(self): # modified function to retain only the filter cols for the H channel
        print('combining modified precursors')
        
        # Define aggregation functions for columns
        keys = self.filter_cols
        agg_functions = {key: 'first' for key in keys}
        agg_functions['intensity'] = 'sum'  # Use 'sum' or 'first' as per your requirement
    
        # Aggregate data using groupby and agg function
        agg_df = self.parsed_df.groupby(['Run', 'Protein.Group', 'Precursor', 'quantity type', 'Label']).agg(agg_functions).reset_index()
    
        # Pivot the table to get 'H', 'M', 'L' intensities in separate columns
        pivoted_df = agg_df.pivot_table(index=['Run', 'Protein.Group', 'Precursor', 'quantity type'], columns='Label', values='intensity', fill_value=0).reset_index()
    
        # Rename columns based on the label
        pivoted_df.columns.name = None
        pivoted_df = pivoted_df.rename(columns={'H': 'H intensity', 'M': 'M intensity', 'L': 'L intensity'})
    
        # Filter the original dataframe for 'H' label and keep only filter columns
        h_filter_df = self.parsed_df[self.parsed_df['Label'] == 'H'][['Run', 'Protein.Group', 'Precursor', 'quantity type'] + self.filter_cols]
    
        # Merge the h_filter_df with the pivoted_df
        df = pd.merge(h_filter_df, pivoted_df, on=['Run', 'Protein.Group', 'Precursor', 'quantity type'])
    
        # Sort and drop duplicates if needed
        df = df.sort_values(['Run', 'Protein.Group', 'Precursor', 'quantity type'], ascending=False).drop_duplicates(subset=['Run', 'Protein.Group', 'Precursor', 'quantity type'])

        return df


    def stack_intensities(self):
        print("Stacking intensities")
        df = self.combined_precursors
        columns_to_check = ['H intensity', 'M intensity', 'L intensity']
        for col in columns_to_check:
            if col not in df.columns:
                df[col] = 0
        df['Precursor.Quantity'] = df['H intensity'] + df['M intensity'] + df['L intensity']
        df['L to stack ratio'] = df['L intensity']/df['Precursor.Quantity']
        df['M to stack ratio'] = df['M intensity']/df['Precursor.Quantity']
        df['H to stack ratio'] = df['H intensity']/df['Precursor.Quantity']
        df = self.drop_nan(df)
        return df


    def drop_nan(self, df):
        cols_to_check = [
        'Precursor.Quantity', 'H intensity', 'L intensity', 'M intensity', 
        'L to stack ratio', 'M to stack ratio', 'H to stack ratio'
        ]
        # Create a mask where all values in the specified columns are either 0 or NaN
        mask = df[cols_to_check].isin([0, np.nan]).all(axis=1)
        
        # Use the mask to drop rows from the DataFrame
        df = df[~mask]
        
        #rename precursor for directLFQ
        df = df.rename(columns={'Precursor':'Precursor.Id'})
        return df[['Run', 'Protein.Group', 'Precursor.Id','Precursor.Quantity','quantity type','H intensity','M intensity','L intensity','H to stack ratio', 'M to stack ratio', 'L to stack ratio']
                  +self.filter_cols]

  