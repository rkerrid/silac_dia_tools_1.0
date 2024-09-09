# -*- coding: utf-8 -*-
"""
Created on Sun Aug 13 13:09:52 2023

@author: rkerrid
"""

## Extracting raw file data 
def extract_values(ms1_df, peptide_details):
    scan_number = peptide_details["scan_number"].values[0]
    scan = ms1_df[ms1_df['scanms1']==scan_number][["mz_values_ms1","intensity_values_ms1", "rtms1"]]

    # Extract lists from DataFrame
    x_values = scan['mz_values_ms1'].tolist()
    y_values = scan['intensity_values_ms1'].tolist()
    x_values = [item for sublist in x_values for item in sublist]
    y_values = [item for sublist in y_values for item in sublist]
   
    return x_values, y_values


##Getting peptide details from MQ files
# Locate the most abundant peptides containing R, K, and P
def get_most_abundant_krp_peptide(msms, evidence, channel):
    msms = msms[msms['Proteins'].notna()]
    df_sorted = msms.sort_values(by='Precursor_Intensity', ascending=False)

    conditions = [
        df_sorted['Sequence'].str.contains('P'),
        (df_sorted['Sequence'].str.contains('K') & ~df_sorted['Sequence'].str.contains('P')),
        (df_sorted['Sequence'].str.contains('R') & ~df_sorted['Sequence'].str.contains('P'))
    ]

    peptides = []
    for condition in conditions:
        sequence = df_sorted[condition].iloc[0]
        df_sorted = df_sorted[df_sorted['Sequence'] != sequence['Sequence']]
        peptides.append(get_peptide_details(sequence, evidence, channel))

    return peptides #list of peptide details


def get_AA_mass(sequence, channel):
    if channel == 'H':
        if "K" in sequence:
            return 8
        elif "R" in sequence:
            return 10
        return 1
    elif channel == 'M':
        if "K" in sequence:
            return 4
        elif "R" in sequence:
            return 6
        return 1

# Locate the most abundant heavy peptide
def get_most_abundant_heavy(msms, evidence, channel):
    msms = msms[msms['Proteins'].notna()]
    msms = msms[msms["Labeling_state"]==1]
    df_sorted = msms.sort_values(by='Precursor_Intensity', ascending=False)
    largest_heavy = df_sorted.iloc[0]
    peptide = get_peptide_details(largest_heavy, evidence, channel)
    return peptide
    

def get_protein_precursors(poi, protein_groups, msms, evidence, channel):
    #check that protein is in df and filter pg groups for peptide ids 
    contains_gene = protein_groups['Protein IDs'].str.contains(poi)
    filtered_proteins = protein_groups[contains_gene]
    
    peptide_ids = filtered_proteins['Peptide IDs']
    peptide_ids = peptide_ids.values.tolist()
    peptide_ids = peptide_ids[0].split(';')
    print('code adjusted')

    # get msms details
    msms_details = ''
    peptides = []
    for peptide_id in peptide_ids:
        msms_details = msms[msms['Peptide_ID']==int(peptide_id)]
        msms_sorted = msms_details.sort_values(by='Precursor_Intensity', ascending=False)
        top_precursor = msms_sorted.iloc[0]
        peptides.append(top_precursor)
    
    peptide_details = []
    # loop through peptides and get peptide details for each one
    for peptide in peptides:
        peptide_details.append(get_peptide_details(peptide, evidence, channel))
    return peptide_details


# Extract peptide details
def get_peptide_details(peptide, evidence, channel):
    full_scan_number = peptide["Precursor_full_scan_number"]
    evidence_id = peptide['Evidence_ID']
    labeling_state = peptide['Labeling_state']
    
    peptide_details = evidence[evidence["id"] == evidence_id][['Raw_file','Sequence', 'm/z', 'Charge', 'Retention_time']]    
    peptide_details['scan_number'] = full_scan_number
    peptide_details["labeling_state"] = labeling_state
    
    sequence = peptide_details["Sequence"].values[0] 
    peptide_details["AA_mass"] = get_AA_mass(sequence, channel)
    
    return peptide_details  







