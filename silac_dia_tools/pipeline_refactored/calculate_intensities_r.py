# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 15:05:24 2024

@author: robbi
"""

def calculate_href_intensities(df):
    df_copy = df.copy(deep=True)
    
    # Calculate median H value and reset index to make it a DataFrame
    h_ref = df_copy.groupby('Protein.Group')['H'].median().reset_index()
    
    # Rename the median column to 'h_ref'
    h_ref = h_ref.rename(columns={'H': 'h_ref'})
    
    # Merge the original DataFrame with the h_ref DataFrame
    merged_df = df.merge(h_ref, on='Protein.Group', how='inner')
    
    # calculate factor to multiply other chanels by dividing href by original H intensity for each PG
    merged_df['factor'] = merged_df['h_ref']/merged_df['H']
    
    # Normalize other chanels with this factor
    merged_df['H_norm'] = merged_df['H']*merged_df['factor']
    merged_df['M_norm'] = merged_df['M']*merged_df['factor']
    merged_df['L_norm'] = merged_df['L']*merged_df['factor']
    
    return merged_df

def output_protein_groups(df, quantification, path):
    if quantification == 'href':
        cols = ['Run', 'Protein.Group', 'H_norm', 'M_norm', 'L_norm']
        df = df[cols]
        df = df.rename(columns={'H_norm': 'H', 'M_norm': 'M', 'L_norm': 'L'})

        # Pivoting for 'H'
        h_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='H')
        
        # Pivoting for 'M'
        m_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='M')
        
        # Pivoting for 'L'
        l_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='L')
        
        # then output each table to csv for h.href, l.href, m.href
        h_pivot_df.to_csv(f'{path}href_href.csv', sep=',')
        m_pivot_df.to_csv(f'{path}nsp_href.csv', sep=',')
        h_pivot_df.to_csv(f'{path}light_href.csv', sep=',')

    return h_pivot_df, m_pivot_df, l_pivot_df

# for dlfq

# def calculate_dlfq_intensities(df):
#     df_copy = df.copy(deep=True)
    
#     # Calculate median H value and reset index to make it a DataFrame
#     h_ref = df_copy.groupby('Protein.Group')['H'].median().reset_index()
    
#     # Rename the median column to 'h_ref'
#     h_ref = h_ref.rename(columns={'H': 'h_ref'})
    
#     # Merge the original DataFrame with the h_ref DataFrame
#     merged_df = df.merge(h_ref, on='Protein.Group', how='inner')
    
#     # calculate factor to multiply other chanels by dividing href by original H intensity for each PG
#     merged_df['factor'] = merged_df['h_ref']/merged_df['H']
    
#     # Normalize other chanels with this factor
#     merged_df['H_norm'] = merged_df['H']*merged_df['factor']
#     merged_df['M_norm'] = merged_df['M']*merged_df['factor']
#     merged_df['L_norm'] = merged_df['L']*merged_df['factor']
    
#     return merged_df

# def output_protein_groups_dlfq(df, quantification, path):
#     if quantification == 'href':
#         cols = ['Run', 'Protein.Group', 'H_norm', 'M_norm', 'L_norm']
#         df = df[cols]
#         df = df.rename(columns={'H_norm': 'H', 'M_norm': 'M', 'L_norm': 'L'})

#         # Pivoting for 'H'
#         h_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='H')
        
#         # Pivoting for 'M'
#         m_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='M')
        
#         # Pivoting for 'L'
#         l_pivot_df = df.pivot(index='Protein.Group', columns='Run', values='L')
        
#         # then output each table to csv for h.href, l.href, m.href
#         h_pivot_df.to_csv(f'{path}href_href.csv', sep=',')
#         m_pivot_df.to_csv(f'{path}nsp_href.csv', sep=',')
#         h_pivot_df.to_csv(f'{path}light_href.csv', sep=',')

#     return h_pivot_df, m_pivot_df, l_pivot_df
