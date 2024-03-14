# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 09:43:38 2024

@author: robbi
"""
### annas R code for tripple silac


  data %>% 
    select(Run, Protein.Group, Precursor.Id, Stripped.Sequence.Charge, Ms1.Translated, Precursor.Translated) %>% 
    mutate(Channel = str_remove(str_extract(Precursor.Id, "-[LMH]"), "-")) %>%  
    select(-Precursor.Id) %>% 
    pivot_longer(contains("Translated"), names_to = "Intensity.Type", values_to = "Intensity") %>% 
    filter(Intensity > 0 & !is.na(Intensity)) %>% 
    pivot_wider(everything(), names_from = Channel, values_from = Intensity) %>% 
    inner_join(filterSet) -> data

  
  if(PGcalculation != "standard"){
    data %>% 
      filter(Intensity.Type == PGcalculation) -> data
  }
  # calculate ratios depending on the channels
  if(numberChannels == 2){
    if("H" %in% colnames(data) & "L" %in% colnames(data)){
      data %>% 
        mutate(LvsH = log10(L/H)) %>% 
        group_by(Run, Protein.Group) %>% 
        # PG level L/H ratio and L abundance normalized to global H
        mutate(LvsH.PG = median(LvsH, na.rm = T)) %>% 
        distinct() -> PG.Ratios
    }else if("H" %in% colnames(data) & "M" %in% colnames(data)){
      data %>% 
        mutate(MvsH = log10(M/H)) %>% 
        group_by(Run, Protein.Group) %>% 
        # PG level L/H ratio and L abundance normalized to global H
        mutate(MvsH.PG = median(MvsH, na.rm = T)) %>% 
        distinct() -> PG.Ratios
    }else{
      data %>% 
        mutate(LvsM = log10(L/M)) %>% 
        group_by(Run, Protein.Group) %>% 
        # PG level L/H ratio and L abundance normalized to global H
        mutate(LvsM.PG = median(LvsM, na.rm = T)) %>% 
        distinct() -> PG.Ratios
    }
   
    

      

  }else if(numberChannels == 3){
    data %>% 
      mutate(LvsH = log10(L/H), MvsH = log10(M/H), LvsM = log10(L/M)) %>% 
      group_by(Run, Protein.Group) %>% 
      # PG level L/H ratio and L abundance normalized to global H
      mutate(LvsH.PG = median(LvsH, na.rm = T), 
             MvsH.PG = median(MvsH, na.rm = T), 
             LvsM.PG = median(LvsM, na.rm = T))  %>% 
      distinct() -> PG.Ratios
  }else{
    stop("Invalid Number of Channels")
  }
    
    
    # add global protein intensity based on reference channel
    PG.Ratios %>% 
      ungroup %>% 
      group_by(Protein.Group) %>% 
      mutate(Global.Log10.Reference.Protein.Intensity = median(log10((!!sym(globalReference))))) %>% 
      ungroup -> PG.Ratios
      
    
    if("NA" %in% colnames(PG.Ratios)){
      PG.Ratios %>% 
        select(-`NA`) -> PG.Ratios
    }
    
    PG.Ratios %>% 
      add_column(use.Filter = useFilter) %>% 
      add_column(precursor.Per.Protein = precursorPerProtein) %>% 
      add_column(global.Reference = globalReference) %>% 
      mutate(across(contains(".PG"), .fns=~.+Global.Log10.Reference.Protein.Intensity, .names = "{.col}.Intensity")) -> PG.Ratios
