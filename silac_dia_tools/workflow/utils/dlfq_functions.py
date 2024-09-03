# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 17:43:03 2023

From the Mann labs https://github.com/MannLabs/directlfq/blob/master/nbdev_nbs/01_lfq_manager.ipynb 
with small edits
"""

#| export
import directlfq.normalization as lfqnorm
import directlfq.protein_intensity_estimation as lfqprot_estimation
import directlfq.utils as lfqutils
import pandas as pd
import directlfq

import warnings


warnings.filterwarnings(action='once')


def run_lfq(input_file, file = '', columns_to_add = [], selected_proteins_file :str = None, mq_protein_groups_txt = None, min_nonan = 1, input_type_to_use = None, maximum_number_of_quadratic_ions_to_use_per_protein = 10, 
number_of_quadratic_samples = 50, num_cores = None, filename_suffix = "", deactivate_normalization = False
):
    """Run the directLFQ pipeline on a given input file. The input file is expected to contain ion intensities. The output is a table containing protein intensities.

    Args:
        input_file (_type_): the input file containing the ion intensities. Usually the output of a search engine.
        columns_to_add (list, optional): additional columns to add to the LFQ intensity output table. They are extraced from the input file. Defaults to [].
        selected_proteins_file (str, optional): if you want to perform normalization only on a subset of proteins, you can pass a .txt file containing the protein IDs, separeted by line breaks. No header expected. Defaults to None.
        mq_protein_groups_txt (_type_, optional): In the case of using MaxQuant data, the proteinGroups.txt table is needed in order to map IDs analogous to MaxQuant. Adding this table improves protein mapping, but is not necessary. Defaults to None.
        min_nonan (int, optional): Min number of ion intensities necessary in order to derive a protein intensity. Increasing the number results in more reliable protein quantification at the cost of losing IDs. Defaults to 1.
        input_type_to_use (_type_, optional): If you want to parse data from the input file in a differing way than specified in the defaults (e.g. extracting MS1 intensities only from a DIANN file), you can name the parsing protocol to be used. The parsing protocols are defined in directlfq/configs/intable_configs.yaml Defaults to None.
        maximum_number_of_quadratic_ions_to_use_per_protein (int, optional): How many ions are used to create the anchor intensity trace (see paper). Increasing might marginally increase performance at the cost of runtime. Defaults to 10.
        number_of_quadratic_samples (int, optional): How many samples are are used to create the anchor intensity trace (see paper). Increasing might marginally increase performance at the cost of runtime. Defaults to 50.
        num_cores (_type_, optional): Num cores to use. Maximum feasible number utilized if set to None. Defaults to None.
    """
    print("Starting directLFQ analysis.")
    input_file = prepare_input_filename(input_file)
    print("reformatting input file, for large files this might take a while.")
    input_file = lfqutils.add_mq_protein_group_ids_if_applicable_and_obtain_annotated_file(input_file, input_type_to_use,mq_protein_groups_txt, columns_to_add)
    input_df = lfqutils.import_data(input_file=input_file, input_type_to_use=input_type_to_use)
    input_df = lfqutils.index_and_log_transform_input_df(input_df)
    input_df = lfqutils.remove_allnan_rows_input_df(input_df)
    
    if not deactivate_normalization:
        print("Performing sample normalization.")
        input_df = lfqnorm.NormalizationManagerSamplesOnSelectedProteins(input_df, num_samples_quadratic=number_of_quadratic_samples, selected_proteins_file=selected_proteins_file).complete_dataframe
    
    print("Estimating lfq intensities.")
    protein_df, ion_df = lfqprot_estimation.estimate_protein_intensities(input_df,min_nonan=min_nonan,num_samples_quadratic=maximum_number_of_quadratic_ions_to_use_per_protein, num_cores = num_cores)
    
    try:
        protein_df = lfqutils.add_columns_to_lfq_results_table(protein_df, input_file, columns_to_add)
    except:
        print("Could not add additional columns to protein table, printing without additional columns.")
    
    print("Writing results files.")
    outfile_basename = get_outfile_basename(input_file, input_type_to_use, selected_proteins_file, deactivate_normalization,filename_suffix)
    save_run_config(outfile_basename, locals())
    save_protein_df(protein_df,file,outfile_basename)
    save_ion_df(ion_df,outfile_basename)
    
    print("Analysis finished!")
    return protein_df

def prepare_input_filename(input_file):
    input_file = fr"{input_file}".replace("\ ", " ").rstrip() #replace escaped spaces with normal spaces and remove trailing whitespace
    return input_file

def get_outfile_basename(input_file, input_type_to_use, selected_proteins_file, deactivate_normalization,filename_suffix):
    outfile_basename = input_file
    outfile_basename += "" if input_type_to_use is None else f".{input_type_to_use}"
    outfile_basename += ".selected_proteins" if selected_proteins_file is not None else ""
    outfile_basename += ".no_norm" if deactivate_normalization else ""
    outfile_basename += filename_suffix
    return outfile_basename

def save_protein_df(protein_df,file, outfile_basename):
    protein_df.to_csv(file, sep = "\t")

def save_ion_df(ion_df, outfile_basename):
    ion_df.to_csv(f"{outfile_basename}.ion_intensities.tsv", sep = "\t")


def save_run_config(outfile_basename, kwargs):
    try:
        df_configs = pd.DataFrame.from_dict(kwargs, orient='index', columns=['value'])
        #add row with directlfq version
        df_configs.loc["directlfq_version"] = directlfq.__version__
        df_configs.to_csv(f"{outfile_basename}.run_config.tsv", sep = "\t")
    except:
        print("Could not save run config.")
        
        
        
        
 