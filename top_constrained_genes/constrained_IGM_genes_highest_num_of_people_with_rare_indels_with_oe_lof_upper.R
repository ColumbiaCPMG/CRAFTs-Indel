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
IGM_indels = fread( file = '20220909_IGM_srAF_sorted.csv', sep = ',', header = TRUE)
IGM_indels$varID = paste0(IGM_indels$CHR, sep="-", IGM_indels$POS, sep="-", IGM_indels$REF,sep="-",IGM_indels$ALT)

ATAV_cols = fread(file = "20220912_IGM_varID_genename_omimdisease_gnomadpli_samplename.csv", sep=",",header = TRUE)
colnames(ATAV_cols) = c("varID", "geneName", "OMIM_disease", "gnomAD_pLI", "sampleName")

## Read in oe_lof_upper column 
oe_lof_upper_varID = fread(file = "20220916_varID_oe_lof_upper_IGM_no_duplic_with_header.csv", sep = ",", header = TRUE)
colnames(oe_lof_upper_varID) = c("varID", "geneName", "oe_lof_upper")
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
constrained_rare_sAF_common_rAF5 = constrained %>% filter(constrained$varID %in% rare_sAF_common_rAF5$varID)
constrained_rare_sAF_common_rAF10 = constrained %>% filter(constrained$varID %in% rare_sAF_common_rAF10$varID)
constrained_rare_sAF_common_rAF15 = constrained %>% filter(constrained$varID %in% rare_sAF_common_rAF15$varID)
constrained_rare_sAF_common_rAF20 = constrained %>% filter(constrained$varID %in% rare_sAF_common_rAF20$varID)

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

## sanity check to make sure that the count is correct to only count the number of distinct sampleNames in each group by gene name 
## success!! matches the dataframe summary_bp5 
test = constrained_rare_sAF_common_rAF5 %>% filter(constrained_rare_sAF_common_rAF5$geneName == "'AAK1'") 
##7
length(unique(test$sampleName))
test2 = constrained_rare_sAF_common_rAF5 %>% filter(constrained_rare_sAF_common_rAF5$geneName == "'ABCA2'") 
##8
length(unique(test2$sampleName))
test3 = constrained_rare_sAF_common_rAF5 %>% filter(constrained_rare_sAF_common_rAF5$geneName == "'ASTN2'") 
##23
length(unique(test3$sampleName))

## sanity check to make sure that the number of genes are retained for bp 5 after removing duplicates 
## the number of unique gene names should be the same 
length(unique(summary_bp5$geneName))
length(unique(OMIM_pLI_bp5$geneName))
length(unique(OMIM_pLI_no_duplic_bp5$geneName))
length(unique(output_5$geneName))
## OMIM_pLI_bp5 dataframe has duplicates so that's why the length is longer
## After removing duplicates in the OMIM_pLI_no_duplic_bp5 dataframe, the length is back to the same as the unique number of gene names 
length((summary_bp5$geneName))
length((OMIM_pLI_bp5$geneName))
length((OMIM_pLI_no_duplic_bp5$geneName))
length((output_5$geneName))

## merge with oe_lof_upper columns
oe_lof_upper_out_5 = left_join(output_5, oe_lof_upper, by="geneName", keep=FALSE)
oe_lof_upper_out_10 = left_join(output_10, oe_lof_upper, by="geneName", keep=FALSE)
oe_lof_upper_out_15 = left_join(output_15, oe_lof_upper, by="geneName", keep=FALSE)
oe_lof_upper_out_20 = left_join(output_20, oe_lof_upper, by="geneName", keep=FALSE)

## write out the tables 
write.table(oe_lof_upper_out_5, "20220926_constrained_highest_num_people_rare_sAF_common_rAF5_oe_lof_upper.tsv",  sep="\t", row.names = FALSE, quote = FALSE)
write.table(oe_lof_upper_out_10, "20220926_constrained_highest_num_people_rare_sAF_common_rAF10_oe_lof_upper.tsv", sep="\t", row.names = FALSE, quote = FALSE)
write.table(oe_lof_upper_out_15, "20220926_constrained_highest_num_people_rare_sAF_common_rAF15_oe_lof_upper.tsv", sep="\t", row.names = FALSE, quote = FALSE)
write.table(oe_lof_upper_out_20, "20220926_constrained_highest_num_people_rare_sAF_common_rAF20_oe_lof_upper.tsv", sep="\t", row.names= FALSE, quote = FALSE)

