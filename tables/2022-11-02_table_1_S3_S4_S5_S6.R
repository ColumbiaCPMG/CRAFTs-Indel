## September 26, 2022
## Sandy Yang
## rAF Project
## This code is to identify constrained (pli > 0.5) OMIM genes with the highest number of people with rare indels based on sAF but are common based on rAF5
## This code includes the oe_lof_upper column  

library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)

rare_IGM = 10^-4 

###### read in files
IGM_indels = fread( file = '2023-03-06_IGM_sAF_rAF_lte_50bp.csv', sep = ',', header = TRUE)

ATAV_cols = fread(file = "20220912_IGM_varID_genename_omimdisease_gnomadpli_samplename.csv", sep=",",header = TRUE)
colnames(ATAV_cols) = c("VarID", "geneName", "OMIM_disease", "gnomAD_pLI", "sampleName")

## Read in oe_lof_upper column 
oe_lof_upper_varID = fread(file = "20220916_varID_oe_lof_upper_IGM_no_duplic_with_header.csv", sep = ",", header = TRUE)
colnames(oe_lof_upper_varID) = c("VarID", "geneName", "oe_lof_upper")
oe_lof_upper = distinct(oe_lof_upper_varID[, c("geneName", "oe_lof_upper")])
######

## get only "constrained indels" or indels on constrained genes 
constrained = ATAV_cols %>% filter (gnomAD_pLI > 0.5)


## get only rare sAF common rAFx indels
rare_sAF_common_rAF5 = IGM_indels %>% filter(IGM_indels$sAF <= rare_IGM & IGM_indels$rAF5 > rare_IGM)
rare_sAF_common_rAF10 = IGM_indels %>% filter(IGM_indels$sAF <= rare_IGM & IGM_indels$rAF10 > rare_IGM)
rare_sAF_common_rAF15 = IGM_indels %>% filter(IGM_indels$sAF <= rare_IGM & IGM_indels$rAF15 > rare_IGM)
rare_sAF_common_rAF20 = IGM_indels %>% filter(IGM_indels$sAF <= rare_IGM & IGM_indels$rAF20 > rare_IGM)

## Get only the constrained indels such that the indels are rare by sAF and common by rAFx 
constrained_rare_sAF_common_rAF5 = constrained %>% filter(constrained$VarID %in% rare_sAF_common_rAF5$VarID)
constrained_rare_sAF_common_rAF10 = constrained %>% filter(constrained$VarID %in% rare_sAF_common_rAF10$VarID)
constrained_rare_sAF_common_rAF15 = constrained %>% filter(constrained$VarID %in% rare_sAF_common_rAF15$VarID)
constrained_rare_sAF_common_rAF20 = constrained %>% filter(constrained$VarID %in% rare_sAF_common_rAF20$VarID)

## summary table with unique gene names with unique sample names as well
summary_bp5 = constrained_rare_sAF_common_rAF5 %>% group_by(geneName) %>% summarise(count=n_distinct(sampleName))
summary_bp10 = constrained_rare_sAF_common_rAF10 %>% group_by(geneName) %>% summarise(count=n_distinct(sampleName))
summary_bp15 = constrained_rare_sAF_common_rAF15 %>% group_by(geneName) %>% summarise(count=n_distinct(sampleName))
summary_bp20 = constrained_rare_sAF_common_rAF20 %>% group_by(geneName) %>% summarise(count=n_distinct(sampleName))

colnames(summary_bp5) = c("geneName", "Number of Unique Individuals with Indels in this Gene")
colnames(summary_bp10) = c("geneName", "Number of Unique Individuals with Indels in this Gene")
colnames(summary_bp15) = c("geneName", "Number of Unique Individuals with Indels in this Gene")
colnames(summary_bp20) = c("geneName", "Number of Unique Individuals with Indels in this Gene")

## Add OMIM and PLI columns to dataframes 
OMIM_pLI_bp5 = constrained_rare_sAF_common_rAF5 %>% select(geneName, OMIM_disease, gnomAD_pLI) %>% group_by(geneName)
OMIM_pLI_bp10 = constrained_rare_sAF_common_rAF10 %>% select(geneName, OMIM_disease, gnomAD_pLI) %>% group_by(geneName)
OMIM_pLI_bp15 = constrained_rare_sAF_common_rAF15 %>% select(geneName, OMIM_disease, gnomAD_pLI) %>% group_by(geneName)
OMIM_pLI_bp20 = constrained_rare_sAF_common_rAF20 %>% select(geneName, OMIM_disease, gnomAD_pLI) %>% group_by(geneName)

## Remove duplicate entries 
OMIM_pLI_no_duplic_bp5 = distinct(OMIM_pLI_bp5)
OMIM_pLI_no_duplic_bp10 = distinct(OMIM_pLI_bp10)
OMIM_pLI_no_duplic_bp15 = distinct(OMIM_pLI_bp15)
OMIM_pLI_no_duplic_bp20 = distinct(OMIM_pLI_bp20)

## merge the summary table with the OMIM PLI columns
output_5 = left_join(summary_bp5, OMIM_pLI_no_duplic_bp5, by = "geneName", keep=FALSE)
output_10 = left_join(summary_bp10, OMIM_pLI_no_duplic_bp10, by = "geneName", keep=FALSE)
output_15 = left_join(summary_bp15, OMIM_pLI_no_duplic_bp15, by = "geneName", keep=FALSE)
output_20 = left_join(summary_bp20, OMIM_pLI_no_duplic_bp20, by= "geneName", keep=FALSE)


## merge with oe_lof_upper columns
oe_lof_upper_out_5 = left_join(output_5, oe_lof_upper, by="geneName", keep=FALSE)
oe_lof_upper_out_10 = left_join(output_10, oe_lof_upper, by="geneName", keep=FALSE)
oe_lof_upper_out_15 = left_join(output_15, oe_lof_upper, by="geneName", keep=FALSE)
oe_lof_upper_out_20 = left_join(output_20, oe_lof_upper, by="geneName", keep=FALSE)




## Table 1
## Genes with > 50 individuals associated with disorders for 10bp sliding window 
## out of these, there are few that are not associated with dominant disorders (6 of these genes)
table1_all = oe_lof_upper_out_5 %>% filter((oe_lof_upper_out_5$`Number of Unique Individuals with Indels in this Gene` > 50) & (!is.na(oe_lof_upper_out_5$OMIM_disease)))
table1_autosomal_dominant_disorders = table1_all[table1_all$OMIM_disease %like% "Autosomal dominant",]

#write.table(table1_all, paste0("/Users/sy3115/Documents/Data/rAF_paper/finalized_project_data/", Sys.Date(), "_table1_all.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

write.table(table1_autosomal_dominant_disorders, paste0( Sys.Date(), "_table1_autosomal_dominant.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## Table S1
table_s1 = oe_lof_upper_out_5 %>% filter(oe_lof_upper_out_5$`Number of Unique Individuals with Indels in this Gene`  > 50)

write.table(table_s1, paste0( Sys.Date(), "_table_s1.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## Table S2
table_s2 = oe_lof_upper_out_10 %>% filter(oe_lof_upper_out_10$`Number of Unique Individuals with Indels in this Gene`  > 50)

write.table(table_s2, paste0( Sys.Date(), "_table_s2.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## Table S3 
table_s3 = oe_lof_upper_out_15 %>% filter(oe_lof_upper_out_15$`Number of Unique Individuals with Indels in this Gene`  > 50)

write.table(table_s3, paste0( Sys.Date(), "_table_s3.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## Table S4 
table_s4 = oe_lof_upper_out_20 %>% filter(oe_lof_upper_out_20$`Number of Unique Individuals with Indels in this Gene`  > 50)

write.table(table_s4, paste0( Sys.Date(), "_table_s4.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## determine if the genes found in the other tables are all in table s4. 
## yes, all the genes in tables 1, s1, s2, s3 are also in s4. 
genes_1_in_s4 = nrow(table1_autosomal_dominant_disorders %>% filter(table1_autosomal_dominant_disorders$geneName %in% table_s4$geneName))
genes_s1_in_s4 = nrow(table_s1 %>% filter(table_s1$geneName %in% table_s4$geneName))
genes_s2_in_s4 = nrow(table_s2 %>% filter(table_s2$geneName %in% table_s4$geneName))
genes_s3_in_s4 = nrow(table_s3 %>% filter(table_s3$geneName %in% table_s4$geneName))


#######################################
####### For the written section #######
#######################################

## XX individuals (XX% of the IGM chort) carry a rare indel in genes associated with autosomal dominant disorders in 10bp window


## merge the atav cols with the IGM_indels 
atav_and_indels = merge(ATAV_cols, IGM_indels, by = "VarID")
## filter for dominant genes (Autosomal dominant)
AD_genes_only = atav_and_indels[atav_and_indels$OMIM_disease %like% "Autosomal dominant",]

## filter for rare_sAF_rare_rAF5 
rare_sAF_common_rAF5_AD = AD_genes_only %>% filter (sAF <= rare_IGM & rAF5 > rare_IGM)
## find unique individuals 
a = unique(rare_sAF_common_rAF5_AD$sampleName)
num_of_ind_with_rare_indel_AD = 8391

b = unique(atav_and_indels$sampleName)
tot_num_samples = 40761

percent_ind_of_IGM_cohort = num_of_ind_with_rare_indel_AD / tot_num_samples * 100 
percent_ind_of_IGM_cohort

## XX individuals carry a suspicious indel in the bp window 
sus_indels_window10 = atav_and_indels %>% filter(sAF <= rare_IGM & rAF5 > rare_IGM)
# find unique individuals 
c = unique(sus_indels_window10$sampleName)
num_of_ind_sus_indels = 28732 

percent_ind_of_IGM_cohort_sus_indels = num_of_ind_sus_indels / tot_num_samples * 100 
percent_ind_of_IGM_cohort_sus_indels
