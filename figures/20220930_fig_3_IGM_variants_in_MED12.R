## Date: September 30, 2022
## Author: Sandy Yang 

## This code is to make a data frame with all the IGM indels in all the windows that includes these this indel in CDKN1C, including their allele count, their impact, and if they are in ClinVar. 
## This code also includes the columns for gene name, effect of the variant and ClinVar classification

## Indel: X-70360679-GGCAGCAGCA-G
## Gene: MED12

## import packages
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)

## define rare 

rare = 10^-4

## define the bp windows
bp_intervals = list("5", "10", "15", "20")

########
## Read in files 
## sAF, rAF and groupID info for IGM indels 
indels = fread( file = '2022-12-30_IGM_sAF_rAF.csv', sep = ',', header = TRUE)
indels$varID = paste0(indels$CHR, sep="-", indels$POS, sep="-", indels$REF,sep="-",indels$ALT)
## Read in ClinVar file 
clinvar = fread(file = 'ClinVar_2022_04_03.tsv', sep = '\t', header = FALSE)
colnames(clinvar) = c("CHR", "POS", "REF", "ALT", "V5", "V6", "V7", "V8", "V9", "predictedDeleteriousness", "V11", "V12", "V13", "V14", "V15", "V16" )
clinvar$varID = paste0(clinvar$CHR, sep="-", clinvar$POS, sep="-", clinvar$REF,sep="-",clinvar$ALT)
## Read in impact information from ATAV
impact = fread(file = "20220916_varID_impact_IGM_no_duplic_with_header.csv", sep = ",", header = TRUE)
colnames(impact) = c("varID", "impact")
## Read in effect file 
effect_gene_name = fread(file="20220916_varID_effect_genename_IGM_no_duplic_no_header.csv", sep = ",", header = TRUE )
colnames(effect_gene_name) = c("varID", "Effect", "geneName")
## Read in oe_lof_upper column 
oe_lof_upper = fread(file = "20220916_varID_oe_lof_upper_IGM_no_duplic_with_header.csv", sep = ",", header = TRUE)
colnames(oe_lof_upper) = c("varID", "geneName", "oe_lof_upper")
## Read in pli column
pli_duplic = fread(file = "20220912_IGM_varID_genename_omimdisease_gnomadpli_samplename.csv", sep=",",header = TRUE)
colnames(pli_duplic) = c("varID", "geneName", "OMIM_disease", "gnomAD_pLI", "sampleName")
pli = distinct(pli_duplic[, c("varID", "gnomAD_pLI")])

########
## identify the indels
indel1 = indels %>% filter(indels$varID == "X-70360679-GGCAGCAGCA-G")
##indel2 = indels %>% filter(indels$varID == "16-2185472-TCA-T")

CHR_indel1 = indel1$CHR
##CHR_indel2 = indel2$CHR

for (bp in bp_intervals) {
  ## Identify the groupIDs and CHR for these two indels
  GROUPID_bp = paste0("GROUPID_bp",bp)
  GROUPID_indel1 = paste0("GROUPID_indel1_bp",bp)
  ##GROUPID_indel2 = paste0("GROUPID_indel2_bp",bp)
  
  assign(GROUPID_indel1, indel1[[get("GROUPID_bp")]])
  ##assign(GROUPID_indel2, indel2[[get("GROUPID_bp")]])
  
  ## Find the indels within these groupIDs and CHR 
  indels_in_GROUPID_indel_1 = paste0("indels_in_GROUPID_", bp, "_indel_1")
  ##indels_in_GROUPID_indel_2 = paste0("indels_in_GROUPID_", bp, "_indel_2")
  
  assign(indels_in_GROUPID_indel_1, indels %>% filter(indels[[get("GROUPID_bp")]] == get(get("GROUPID_indel1")) & indels$CHR == CHR_indel1))
  ##assign(indels_in_GROUPID_indel_2, indels %>% filter(indels[[get("GROUPID_bp")]] == get(get("GROUPID_indel2")) & indels$CHR == CHR_indel2))

  
  ## Get a data frame with allele count, impact, whether they are in ClinVar 
  ## merge with Impact column
  
  impact_indel1 = paste0("impact_indels_in_GROUPID_bp", bp, "_indel1")
  ##impact_indel2 = paste0("impact_indels_in_GROUPID_bp", bp, "_indel2")

  assign(impact_indel1, merge(get(indels_in_GROUPID_indel_1),impact, by = "varID"))
  ##assign(impact_indel2, merge(get(indels_in_GROUPID_indel_2),impact, by = "varID"))

  ## merge with the effect and gene name 
  effect_indel1 =  paste0("effect_gene_impact_indels_in_GROUPID_bp", bp, "_indel1") 
  ##effect_indel2 =  paste0("effect_gene_impact_indels_in_GROUPID_bp", bp, "_indel2")
  
  assign(effect_indel1, merge( get(impact_indel1), effect_gene_name, by = "varID"))
  ##assign(effect_indel2, merge( get(impact_indel2), effect_gene_name, by = 'varID'))

  ## merge with oe_lof_upper column
  oe_lof_upper_indel1 = paste0("oe_lof_upper_effect_gene_impact_indel1_GROUPID_bp", bp)
  ##oe_lof_upper_indel2 = paste0("oe_lof_upper_effect_gene_impact_indel2_GROUPID_bp", bp)
    
  assign( oe_lof_upper_indel1, merge(get(effect_indel1), oe_lof_upper, by = "varID"))
  ##assign( oe_lof_upper_indel2, merge(get(effect_indel2), oe_lof_upper, by = "varID"))
  
  ## merge with Clinvar classification
  clinvar_indel1 = paste0("clinvar_effect_gene_impact_indels_in_GROUPID_bp",bp, "_indel1")
  ##clinvar_indel2 = paste0("clinvar_effect_gene_impact_indels_in_GROUPID_bp",bp, "_indel2")
  
  assign(clinvar_indel1, left_join(get(oe_lof_upper_indel1), clinvar, by = "varID", keep = FALSE))
  ##assign(clinvar_indel2, left_join(get(oe_lof_upper_indel2), clinvar, by = "varID", keep = FALSE))
  
  ## merge with pLI classification
  pli_indel1 = paste0("pli_clinvar_effect_gene_impact_indels_in_GROUPID_bp",bp,"_indel1")
  ##pli_indel2 = paste0("pli_clinvar_effect_gene_impact_indels_in_GROUPID_bp",bp,"_indel2")
  
  assign(pli_indel1, left_join(get(clinvar_indel1), pli, by ="varID", keep = FALSE))
  ##assign(pli_indel2, left_join(get(clinvar_indel2), pli, by = "varID", keep = FALSE))
  
  ## determine if the indel is in ClinVar and create output 
  output_indel1 = paste0("output_indel1_GROUPID_bp", bp)
  ##output_indel2 = paste0("output_indel2_GROUPID_bp", bp)
  
  assign(output_indel1 ,
    get(pli_indel1) %>% mutate(
      varID = varID,
      Impact = impact,
      Effect = Effect,
      pLI = gnomAD_pLI,
      in_ClinVar = case_when(get(pli_indel1)$varID %in% clinvar$varID ~ "TRUE", TRUE ~ "FALSE"),
      ClinVar_classification = predictedDeleteriousness,
      gene_name = geneName.x,
      oe_lof_upper = oe_lof_upper,
      AC = AC,
      GROUPID_bp5 = GROUPID_bp5,
      GROUPID_bp10 = GROUPID_bp10,
      GROUPID_bp15 = GROUPID_bp15,
      GROUPID_bp20 = GROUPID_bp20,
      AC_bp5 = AlleleCount_bp5,
      AC_bp10 = AlleleCount_bp10,
      AC_bp15 = AlleleCount_bp15,
      AC_bp20 = AlleleCount_bp20,
      sAF = sAF,
      AN = AN, 
      rAF5 = rAF5,
      rAF10 = rAF10, 
      rAF15 = rAF15,
      rAF20 = rAF20,
      CHR = CHR.x,
      POS = POS.x,
      REF = REF.x,
      ALT = ALT.x
    ) %>% select(varID, Impact, Effect, pLI, in_ClinVar, ClinVar_classification, oe_lof_upper, gene_name, AN, AC, AC_bp5, AC_bp10, AC_bp15, AC_bp20, sAF, rAF5, rAF10, rAF15, rAF20, GROUPID_bp5, GROUPID_bp10, GROUPID_bp15, GROUPID_bp20, CHR, POS, REF, ALT) )
  
  ## export data frames
  write.table(get(output_indel1), paste0("20220930_ALL_IGM_indels_in_groupID_bp", bp, "_with_X-70360679-GGCAGCAGCA-G.csv"), sep = ",", row.names = FALSE, quote = FALSE)
  
}
  
  
  

