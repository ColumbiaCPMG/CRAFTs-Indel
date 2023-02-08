### Sandy Yang
### December 8, 2022 
### rAF/Indels Projet
### This code will be to find the numbers for the ClinVar classification analysis of the written results section for gnomAD 
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
## Read in all the files 
IGM_sAF_rAF = fread("2022-11-22_IGM_rsaf_zeroremoved.csv", sep = ",", header = TRUE)
gnomad_sAF_rAF = fread("20221014_gnomad_rsaf_zeroremoved.csv", sep = ",", header = TRUE)
clinvar = fread('ClinVar_2022_04_03.tsv', sep='\t', header = FALSE)

## Set rareness and categorizations 
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
IGM_sAF_rAF$VarID = paste0(IGM_sAF_rAF$CHR, "\t", IGM_sAF_rAF$POS, "\t", IGM_sAF_rAF$REF, "\t", IGM_sAF_rAF$ALT)
clinvar$VarID = paste0(clinvar$V1, "\t", clinvar$V2, "\t", clinvar$V3, "\t", clinvar$V4)

## Merge sAF and rAF with clinvar and ONLY KEEP THE ROWS THAT ARE IN CLINVAR 
gnomad_in_clinvar = merge(gnomad_sAF_rAF, clinvar, by = "VarID") ## 27559
IGM_in_clinvar = merge(IGM_sAF_rAF, clinvar, by = "VarID") ## 13217 

## Remove any duplicates, if any 
gnomad_in_clinvar = distinct(gnomad_in_clinvar$VarID) ## 27559
IGM_in_clinvar = distinct(IGM_in_clinvar$VarID) #13217

###########################################################################
## How many rare (sAF <= 10^-4) indels are benign, likely benign, patho and likely patho?
###########################################################################
## gnomAD 
gnomad_rare_sAF = gnomad_in_clinvar %>% filter(sAF <= rare)

gnomad_rare_sAF_benign = nrow(gnomad_rare_sAF %>% filter(V10 %in% benign))
gnomad_rare_sAF_LB = nrow(gnomad_rare_sAF %>% filter(V10 %in% LB))
gnomad_rare_sAF_patho = nrow(gnomad_rare_sAF %>% filter(V10 %in% patho))
gnomad_rare_sAF_LP = nrow(gnomad_rare_sAF %>% filter(V10 %in% LP))

## IGM
IGM_rare_sAF = IGM_in_clinvar %>% filter(sAF <= rare)

IGM_rare_sAF_benign = nrow(IGM_rare_sAF %>% filter(V10 %in% benign))
IGM_rare_sAF_LB = nrow(IGM_rare_sAF %>% filter(V10 %in% LB))
IGM_rare_sAF_patho = nrow(IGM_rare_sAF %>% filter(V10 %in% patho))
IGM_rare_sAF_LP = nrow(IGM_rare_sAF %>% filter(V10 %in% LP))

###########################################################################
## To calculate percentage of LP_patho and LB_benign that are common by rAF5 (10bp window) for gnomAD
###########################################################################
## gnomAD
gnomad_patho = gnomad_in_clinvar %>% filter(V10 %in% patho)
gnomad_LP = gnomad_in_clinvar %>% filter(V10 %in% LP)
gnomad_benign = gnomad_in_clinvar %>% filter(V10 %in% benign)
gnomad_LB = gnomad_in_clinvar %>% filter(V10 %in% LB)

gnomad_patho_LP_com_rAF5 = (nrow(gnomad_patho %>% filter(rAF5 > rare)) + nrow(gnomad_LP %>% filter(rAF5 > rare))) / (nrow(gnomad_patho) + nrow(gnomad_LP))
gnomad_benign_LB_com_rAF5 = (nrow(gnomad_benign %>% filter(rAF5 > rare)) + nrow(gnomad_LB %>% filter(rAF5 > rare))) / (nrow(gnomad_benign) + nrow(gnomad_LB))

gnomad_patho_LP_com_rAF5
gnomad_benign_LB_com_rAF5

###########################################################################
## To calculate percentage of LP_patho and LB_benign that are common by rAF5 (10bp window) for IGM
###########################################################################
## IGM
IGM_patho = IGM_in_clinvar %>% filter(V10 %in% patho)
IGM_LP = IGM_in_clinvar %>% filter(V10 %in% LP)
IGM_benign = IGM_in_clinvar %>% filter(V10 %in% benign)
IGM_LB = IGM_in_clinvar %>% filter(V10 %in% LB)

IGM_patho_LP_com_rAF5 = (nrow(IGM_patho %>% filter(rAF5 > rare)) + nrow(IGM_LP %>% filter(rAF5 > rare))) / (nrow(IGM_patho) + nrow(IGM_LP))
IGM_benign_LB_com_rAF5 = (nrow(IGM_benign %>% filter(rAF5 > rare)) + nrow(IGM_LB %>% filter(rAF5 > rare))) / (nrow(IGM_benign) + nrow(IGM_LB))

IGM_patho_LP_com_rAF5
IGM_benign_LB_com_rAF5

###########################################################################
## Find num of rare sAF indels classified as patho/LP that has rAF>1% 
###########################################################################
## gnomAD
patho_LP_rare_gnomad_rAF_gt1 = gnomad_rare_sAF %>% filter(V10 %in% c(patho, LP)) %>% filter(rAF5 > 0.01)
nrow(patho_LP_rare_gnomad_rAF_gt1) 
## IGM
patho_LP_rare_IGM_rAF_gt1 = IGM_rare_sAF %>% filter(V10 %in% c(patho, LP)) %>% filter(rAF5 > 0.01)
nrow(patho_LP_rare_IGM_rAF_gt1)





