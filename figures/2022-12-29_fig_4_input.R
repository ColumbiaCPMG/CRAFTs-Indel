### Sandy Yang
### December 29, 2022 
### rAF/Indels Projet
### This code will be to create figure 4 csv input 
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)

## Read in all the files 
IGM_sAF_rAF = fread("2023-03-06_IGM_sAF_rAF_lte_50bp.csv", sep = ",", header = TRUE)
gnomad_sAF_rAF = fread("2023-03-06_gnomAD_sAF_rAF_lte_50bp.csv", sep = ",", header = TRUE)
clinvar = fread('ClinVar_2022_04_03.tsv', sep='\t', header = FALSE)


## Set rareness, rAF5_graph_threshold,  and categorizations 
rare = 10^-4

benign = c('Benign', 'Benign/Likely_benign', 'Benign|_drug_response')
LB = c('Likely_benign', 'Likely_benign|_drug_response|_other',
      'Likely_benign|_other')
conflict = c('Conflicting_interpretations_of_pathogenicity',
            'Conflicting_interpretations_of_pathogenicity|_association',
            'Conflicting_interpretations_of_pathogenicity|_drug_response',
            'Conflicting_interpretations_of_pathogenicity|_drug_response|_other',
            'Conflicting_interpretations_of_pathogenicity|_other',
            'Conflicting_interpretations_of_pathogenicity|_risk_factor',
            'Uncertain_significance',
            'Uncertain_significance|_risk_factor')  #conflict and uncertain significance
LP = c('Likely_pathogenic',
      'Likely_pathogenic|_risk_factor')
patho = c('Pathogenic',
         'Pathogenic/Likely_pathogenic',
         'Pathogenic/Likely_pathogenic|_other', 'Pathogenic|_Affects',
         'Pathogenic|_other', 'Pathogenic|_protective',
         'Pathogenic|_risk_factor')
other = c('Affects', '\\N', 'association',
         'confers_sensitivity', 'drug_response', 'not_provided', 'other',
         'protective', 'risk_factor')

## Give all the files a VarID col 
IGM_sAF_rAF$VarID = paste0(IGM_sAF_rAF$CHR, "-", IGM_sAF_rAF$POS, "-", IGM_sAF_rAF$REF, "-", IGM_sAF_rAF$ALT)
clinvar$VarID = paste0(clinvar$V1, "-", clinvar$V2, "-", clinvar$V3, "-", clinvar$V4)

## Merge sAF and rAF with clinvar and ONLY KEEP THE ROWS THAT ARE IN CLINVAR 
gnomad_in_clinvar = merge(gnomad_sAF_rAF, clinvar, by = "VarID") ## 27559
IGM_in_clinvar = merge(IGM_sAF_rAF, clinvar, by = "VarID") ## 13217 

## Remove any duplicates, if any 
gnomad_in_clinvar = distinct(gnomad_in_clinvar) ## 27559
IGM_in_clinvar = distinct(IGM_in_clinvar) #13217

# sanity check to make sure there is the same number of UNIQUE VarIDs
uniq_varID_gnomad = unique(gnomad_in_clinvar$VarID) #27559 
uniq_varID_clinvar = unique(IGM_in_clinvar$VarID) #13217

## To create the bargraph, we are interested in only looking at rAF5 and RARE SAF indels 
## We need to know how many are benign/likely benign and pathogenic/likely pathogenic  
## Our X-axis is the rAF5 
## Our Y-axis is the number of indels on a log scale of that rAF value 
## Two different bars represent whether or not these are benign/likely benign and pathogenic/likely pathogenic 
## gnomAD 
gnomad_B_LB = gnomad_in_clinvar %>% filter(V10 %in% benign | V10 %in% LB) %>% filter (sAF <= rare)
gnomad_P_LP = gnomad_in_clinvar %>% filter(V10 %in% patho | V10 %in% LP) %>% filter (sAF <= rare) 


## IGM
IGM_B_LB = IGM_in_clinvar %>% filter(V10 %in% benign | V10 %in% LB) %>% filter (sAF <= rare) 
IGM_P_LP = IGM_in_clinvar %>% filter(V10 %in% patho | V10 %in% LP) %>% filter (sAF <= rare)

## Give benign/LB a B/LB label
## Give patho/LP a P/LP label 

#gnomAD
gnomad_B_LB$clinvar_class = "B/LB"
gnomad_P_LP$clinvar_class = "P/LP"

#IGM
IGM_B_LB$clinvar_class = "B/LB"
IGM_P_LP$clinvar_class = "P/LP"

## Give benign/LB a 0 log_key
## Give patho/LP a 1 log_key 
gnomad_B_LB$log_key = 0 #2273
gnomad_P_LP$log_key = 1 #7831

#IGM
IGM_B_LB$log_key = 0 #447
IGM_P_LP$log_key = 1 #4435


## Since we want to make a graph, we need our dataframe to have three columns: clinvar_classification, raf5, num of indels on a log scale with same rAF and clinvar classification
gnomad_output = rbind(gnomad_B_LB, gnomad_P_LP) %>% mutate (
  clinvar_class = clinvar_class,
  rAF5 = rAF5, 
  log_key = log_key
) %>% select (clinvar_class, rAF5, log_key)

IGM_output = rbind(IGM_B_LB, IGM_P_LP) %>% mutate (
  clinvar_class = clinvar_class,
  rAF5 = rAF5,
  log_key = log_key
) %>% select (clinvar_class, rAF5, log_key)

## for these checks, make sure to set gnomad_output and IGM_output to before mutating the column
uniq_varID_gnomad_2 = unique(gnomad_output$VarID) #9258
uniq_varID_igm_2 = unique(IGM_output$VarID) #4339

## output the csv files 

fwrite(gnomad_output, paste0( Sys.Date(), "_figure4_gnomad_input.csv"))
fwrite(IGM_output, paste0( Sys.Date(), "_figure4_IGM_input.csv"))

################################




