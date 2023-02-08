### Sandy Yang
### December 8, 2022 
### rAF/Indels Projet
### This code will be to find the numbers for Supplmentary Table 9 
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
## How many indels in different classifications have a common rAF for each bp window?
## Table S9 
###########################################################################
## gnomAD classifications in ClinVar
gnomad_patho = gnomad_in_clinvar %>% filter(V10 %in% patho)
gnomad_LP = gnomad_in_clinvar %>% filter(V10 %in% LP)
gnomad_benign = gnomad_in_clinvar %>% filter(V10 %in% benign)
gnomad_LB = gnomad_in_clinvar %>% filter(V10 %in% LB)

## window 10 
gnomad_patho_com_rAF5 =nrow(gnomad_patho %>% filter(rAF5 > rare))
gnomad_LP_com_rAF5 =nrow(gnomad_LP %>% filter(rAF5 > rare))
gnomad_benign_com_rAF5 =nrow(gnomad_benign %>% filter(rAF5 > rare))
gnomad_LB_com_rAF5 =nrow(gnomad_LB %>% filter(rAF5 > rare))

gnomad_patho_com_rAF5
gnomad_LP_com_rAF5
gnomad_benign_com_rAF5
gnomad_LB_com_rAF5

## window 20
gnomad_patho_com_rAF10 =nrow(gnomad_patho %>% filter(rAF10 > rare))
gnomad_LP_com_rAF10 =nrow(gnomad_LP %>% filter(rAF10 > rare))
gnomad_benign_com_rAF10 =nrow(gnomad_benign %>% filter(rAF10 > rare))
gnomad_LB_com_rAF10 =nrow(gnomad_LB %>% filter(rAF10 > rare))

gnomad_patho_com_rAF10
gnomad_LP_com_rAF10
gnomad_benign_com_rAF10
gnomad_LB_com_rAF10

## window 30
gnomad_patho_com_rAF15 =nrow(gnomad_patho %>% filter(rAF15 > rare))
gnomad_LP_com_rAF15 =nrow(gnomad_LP %>% filter(rAF15 > rare))
gnomad_benign_com_rAF15 =nrow(gnomad_benign %>% filter(rAF15 > rare))
gnomad_LB_com_rAF15 =nrow(gnomad_LB %>% filter(rAF15 > rare))

gnomad_patho_com_rAF15
gnomad_LP_com_rAF15
gnomad_benign_com_rAF15
gnomad_LB_com_rAF15

## window 40
gnomad_patho_com_rAF20 =nrow(gnomad_patho %>% filter(rAF20 > rare))
gnomad_LP_com_rAF20 =nrow(gnomad_LP %>% filter(rAF20 > rare))
gnomad_benign_com_rAF20 =nrow(gnomad_benign %>% filter(rAF20 > rare))
gnomad_LB_com_rAF20 =nrow(gnomad_LB %>% filter(rAF20 > rare))

gnomad_patho_com_rAF20
gnomad_LP_com_rAF20
gnomad_benign_com_rAF20
gnomad_LB_com_rAF20

###########################################################################
## How many indels in different classifications that have a common rAF?
## Table S9 
###########################################################################
## IGM classifications in ClinVar 
IGM_patho = IGM_in_clinvar %>% filter(V10 %in% patho)
IGM_LP = IGM_in_clinvar %>% filter(V10 %in% LP)
IGM_benign = IGM_in_clinvar %>% filter(V10 %in% benign)
IGM_LB = IGM_in_clinvar %>% filter(V10 %in% LB)

## window 10
IGM_patho_com_rAF5 = nrow(IGM_patho %>% filter(rAF5 > rare))
IGM_LP_com_rAF5 = nrow(IGM_LP %>% filter(rAF5 > rare)) 
IGM_benign_com_rAF5 = nrow(IGM_benign %>% filter( rAF5 > rare))
IGM_LB_com_rAF5 = nrow(IGM_LB %>% filter( rAF5 > rare))

IGM_patho_com_rAF5
IGM_LP_com_rAF5
IGM_benign_com_rAF5
IGM_LB_com_rAF5

## window 20
IGM_patho_com_rAF10 = nrow(IGM_patho %>% filter(rAF10 > rare))
IGM_LP_com_rAF10 = nrow(IGM_LP %>% filter(rAF10 > rare)) 
IGM_benign_com_rAF10 = nrow(IGM_benign %>% filter( rAF10 > rare))
IGM_LB_com_rAF10 = nrow(IGM_LB %>% filter( rAF10 > rare))

IGM_patho_com_rAF10
IGM_LP_com_rAF10
IGM_benign_com_rAF10
IGM_LB_com_rAF10

## window 30
IGM_patho_com_rAF15 = nrow(IGM_patho %>% filter(rAF15 > rare))
IGM_LP_com_rAF15 = nrow(IGM_LP %>% filter(rAF15 > rare)) 
IGM_benign_com_rAF15 = nrow(IGM_benign %>% filter( rAF15 > rare))
IGM_LB_com_rAF15 = nrow(IGM_LB %>% filter( rAF15 > rare))

IGM_patho_com_rAF15
IGM_LP_com_rAF15
IGM_benign_com_rAF15
IGM_LB_com_rAF15

## window 40
IGM_patho_com_rAF20 = nrow(IGM_patho %>% filter(rAF20 > rare))
IGM_LP_com_rAF20 = nrow(IGM_LP %>% filter(rAF20 > rare)) 
IGM_benign_com_rAF20 = nrow(IGM_benign %>% filter( rAF20 > rare))
IGM_LB_com_rAF20 = nrow(IGM_LB %>% filter( rAF20 > rare))

IGM_patho_com_rAF20
IGM_LP_com_rAF20
IGM_benign_com_rAF20
IGM_LB_com_rAF20





