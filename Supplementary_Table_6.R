library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)

setwd("/Users/sy3115/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/rAF_scripts_and_inputs/PublicationReady_Data/")

## Read in file

#Define your threshold:
threshold = (1 * 10^-4)

## define functionsmed_passed_conditions = function(ContaminationPercentage, MeanCoverage, CCDSBasesCov10X, SelfDeclGender, seqGender, XYAvgCovRatio, sample_status) {
# bins with rAF hi P/LP (ONLY 1 RAFHI P/LP) and others
PLP_others_1 = function(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) {
  return (
    (P_rAFhi + LP_rAFhi == 1) & (LB > 0 | B > 0 | VUS > 0) 
  )
}

# bins with rAF hi P/LP (MORE THAN 1 RAFHI P/LP) and others 
PLP_others_gt1 = function(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) {
  return (
    (P_rAFhi + LP_rAFhi > 1) & (LB > 0 | B > 0 | VUS > 0)
  )
}

# bins with P/LP (ONLY 1) only
PLP_only_1 = function(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) {
    (P_rAFhi + LP_rAFhi == 1) & (LB == 0 & B == 0 & VUS == 0)
}

# bins with P/LP (MORE THAN 1) only
PLP_only_gt1 = function(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) {
  return (
    (P_rAFhi + LP_rAFhi > 1) & (LB == 0 & B == 0 & VUS == 0)
  )
}


gnomAD_clinvar = fread("gnomAD/2024-03-23_gnomAD_varIDs_clinvar.csv") %>% select(GID_bp10, sAF, rAF_bp10, VarID, V10)
IGM_clinvar = fread("IGM/2024-03-23_IGM_varIDs_clinvar.csv")  %>% select(GID_bp10, sAF, rAF_bp10, VarID, V10)
UKBB_clinvar = fread("UKBB/2024-03-23_UKBB_varIDs_clinvar.csv")  %>% select(GID_bp10, sAF, rAF_bp10, VarID, V10)

## Define classifications 
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


gnomAD_rAFhi = gnomAD_clinvar %>% filter(sAF <= threshold & rAF_bp10 > threshold)
IGM_rAFhi = IGM_clinvar %>% filter(sAF <= threshold & rAF_bp10 > threshold)
UKBB_rAFhi = UKBB_clinvar %>% filter(sAF <= threshold & rAF_bp10 > threshold)

###################################################################################################
###################################################################################################
###################################################################################################

### now make a df that is bin ID | # P | # LP | # VUS | # LB | # B | # LP_rAFhi | # P_rAFhi 

## label each variant as P/LP/VUS/LB/B/LP_rAFhi/P_rAFhi
gnomAD_table_rAFhi = gnomAD_clinvar %>% mutate (Classification = case_when(V10 %in% benign ~ "B", V10 %in% LB ~ "LB", V10 %in% conflict ~ "VUS", ((V10 %in% LP) & !(VarID %in% gnomAD_rAFhi$VarID)) ~ "LP", ((V10 %in% patho) & !(VarID %in% gnomAD_rAFhi$VarID)) ~ "P", ((V10 %in% LP) & (VarID %in% gnomAD_rAFhi$VarID)) ~ "LP_rAFhi", ((V10 %in% patho) & (VarID %in% gnomAD_rAFhi$VarID)) ~ "P_rAFhi"))

IGM_table_rAFhi = IGM_clinvar %>% mutate (Classification = case_when(V10 %in% benign ~ "B", V10 %in% LB ~ "LB", V10 %in% conflict ~ "VUS", ((V10 %in% LP) & !(VarID %in% IGM_rAFhi$VarID)) ~ "LP", ((V10 %in% patho) & !(VarID %in% IGM_rAFhi$VarID)) ~ "P", ((V10 %in% LP) & (VarID %in% IGM_rAFhi$VarID)) ~ "LP_rAFhi", ((V10 %in% patho) & (VarID %in% IGM_rAFhi$VarID)) ~ "P_rAFhi"))

UKBB_table_rAFhi = UKBB_clinvar %>% mutate (Classification = case_when(V10 %in% benign ~ "B", V10 %in% LB ~ "LB", V10 %in% conflict ~ "VUS", ((V10 %in% LP) & !(VarID %in% UKBB_rAFhi$VarID)) ~ "LP", ((V10 %in% patho) & !(VarID %in% UKBB_rAFhi$VarID)) ~ "P", ((V10 %in% LP) & (VarID %in% UKBB_rAFhi$VarID)) ~ "LP_rAFhi", ((V10 %in% patho) & (VarID %in% UKBB_rAFhi$VarID)) ~ "P_rAFhi"))

### matrix 

gnomAD_GID_matrix_rAFhi = as.data.frame.matrix(table(gnomAD_table_rAFhi$GID_bp10, gnomAD_table_rAFhi$Classification))
IGM_GID_matrix_rAFhi = as.data.frame.matrix(table(IGM_table_rAFhi$GID_bp10, IGM_table_rAFhi$Classification))
UKBB_GID_matrix_rAFhi = as.data.frame.matrix(table(UKBB_table_rAFhi$GID_bp10, UKBB_table_rAFhi$Classification))

## Find the total number of bins with at least one P/LP rAF hi indel 
gnomAD_tot_rAFhi_bins_PLP = gnomAD_GID_matrix_rAFhi %>% filter(P_rAFhi > 0 | LP_rAFhi > 0)
IGM_tot_rAFhi_bins_PLP = IGM_GID_matrix_rAFhi %>% filter(P_rAFhi > 0 | LP_rAFhi > 0)
UKBB_tot_rAFhi_bins_PLP = UKBB_GID_matrix_rAFhi %>% filter( P_rAFhi > 0 | LP_rAFhi > 0)

# rAF bins with rAF hi P/LP (ONLY 1 RAFHI P/LP) and others
gnomAD_bins_PLP_rAFhi_others_1 = gnomAD_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_others_1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)
IGM_bins_PLP_rAFhi_others_1 = IGM_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_others_1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)
UKBB_bins_PLP_rAFhi_others_1 = UKBB_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_others_1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)

# rAF bins with rAF hi P/LP (MORE THAN 1 RAFHI P/LP) and others 
gnomAD_bins_PLP_rAFhi_others_GT1 = gnomAD_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_others_gt1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)
IGM_bins_PLP_rAFhi_others_GT1 = IGM_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_others_gt1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)
UKBB_bins_PLP_rAFhi_others_GT1 = UKBB_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_others_gt1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)

# rAF bins with P/LP (ONLY 1) only
gnomAD_bins_PLP_rAFhi_only_1 = gnomAD_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_only_1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)
IGM_bins_PLP_rAFhi_only_1 = IGM_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_only_1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)
UKBB_bins_PLP_rAFhi_only_1 = UKBB_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_only_1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)

# rAF bins with P/LP (MORE THAN 1) only
gnomAD_bins_PLP_rAFhi_only_GT1 = gnomAD_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_only_gt1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)
IGM_bins_PLP_rAFhi_only_GT1 = IGM_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_only_gt1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)
UKBB_bins_PLP_rAFhi_only_GT1 = UKBB_tot_rAFhi_bins_PLP %>% mutate(keep = case_when(PLP_only_gt1(P_rAFhi, LP_rAFhi, B, LB, VUS, P, LP) ~ T, T ~ F)) %>% filter(keep == T)

#### CREATE A TABLE 

df = c("gnomAD", "IGM", "UKBB")
tot = c(nrow(gnomAD_tot_rAFhi_bins_PLP), nrow(IGM_tot_rAFhi_bins_PLP), nrow(UKBB_tot_rAFhi_bins_PLP))
PLP_1_w_others = c(nrow(gnomAD_bins_PLP_rAFhi_others_1), nrow(IGM_bins_PLP_rAFhi_others_1), nrow(UKBB_bins_PLP_rAFhi_others_1))
PLP_gt1_w_others = c(nrow(gnomAD_bins_PLP_rAFhi_others_GT1), nrow(IGM_bins_PLP_rAFhi_others_GT1), nrow(UKBB_bins_PLP_rAFhi_others_GT1))
PLP_1_only = c(nrow(gnomAD_bins_PLP_rAFhi_only_1), nrow(IGM_bins_PLP_rAFhi_only_1), nrow(UKBB_bins_PLP_rAFhi_only_1))
PLP_gt1_only = c(nrow(gnomAD_bins_PLP_rAFhi_only_GT1), nrow(IGM_bins_PLP_rAFhi_only_GT1), nrow(UKBB_bins_PLP_rAFhi_only_GT1))

table_out = data.frame(df, tot, PLP_1_w_others, PLP_gt1_w_others, PLP_1_only, PLP_gt1_only)

table_out$prct_PLP_1_w_others = paste0(round(PLP_1_w_others/tot, 2) * 100 , "%")
table_out$prct_PLP_gt1_w_others = paste0(round(PLP_gt1_w_others/tot, 2) * 100, "%")
table_out$prct_PLP_1_only = paste0(round(PLP_1_only/tot, 2) * 100, "%")
table_out$prct_PLP_gt1_only = paste0(round(PLP_gt1_only/tot, 2) * 100, "%")

#a = gnomAD_tot_rAFhi_bins_PLP %>% filter(!(row.names(gnomAD_tot_rAFhi_bins_PLP) %in% row.names(gnomAD_bins_PLP_rAFhi_others_1)) & !(row.names(gnomAD_tot_rAFhi_bins_PLP) %in% row.names(gnomAD_bins_PLP_rAFhi_others_GT1)) & !(row.names(gnomAD_tot_rAFhi_bins_PLP) %in% row.names(gnomAD_bins_PLP_rAFhi_only_1)) & !(row.names(gnomAD_tot_rAFhi_bins_PLP) %in% row.names(gnomAD_bins_PLP_rAFhi_only_GT1))) 

output = table_out %>% select(df, tot, PLP_1_w_others, prct_PLP_1_w_others, PLP_gt1_w_others, prct_PLP_gt1_w_others, PLP_1_only, prct_PLP_1_only, PLP_gt1_only, prct_PLP_gt1_only )
colnames(output) = c("Database", "Total bins with at least 1 rAF-hi P/LP indels", "Bins with 1 rAF-hi P/LP indel and other indels", "Percent of bins with 1 rAF-hi P/LP indel and other indels", "Bins with more than 1 rAF-hi P/LP indels and other indels", "Percent of bins with more than 1 rAF-hi P/LP indels and other indels", "Bins with only 1 rAF-hi P/LP indel", "Percent of bins with only 1 rAF-hi P/LP indel", "Bins with more than 1 rAF-hi P/LP indels and no other types of indels", "Percent of bins with more than 1 rAF-hi P/LP indels and no other types of indels")

fwrite(output, "/Users/sy3115/Documents/Data/rAF_paper/finalized_project_data/PublicationReady/percentage_of_reclassification.csv", row.names = F, quote = F)





