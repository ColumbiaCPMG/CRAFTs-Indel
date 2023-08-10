library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)

rare = 10^-4
bp_range = c("10", "20", "30", "40")

input_df1 = fread("/Users/sy3115/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/rAF_scripts_and_inputs/PublicationReady_Data/gnomAD/gnomad.exomes.r2.1.1.sites_indelsonly_rAF_lt50bp.csv", header = TRUE)
input_df2 = fread("/Users/sy3115/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/rAF_scripts_and_inputs/PublicationReady_Data/IGM/2023-03-23_IGM_n39367_indelsonly_rAF_lt50bp.csv", header = TRUE)

clinvar = fread('/Users/sy3115/Desktop/ClinVar_2023_03_18.tsv', sep='\t', header = FALSE, quote="")
clinvar$VarID = paste0(clinvar$V1, "-", clinvar$V2, "-", clinvar$V3, "-", clinvar$V4)
clinvar = clinvar %>% select("VarID", "V10")

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


effects = fread("/Users/sy3115/Documents/Data/rAF_paper/finalized_project_data/final_dev/2023-04-21_IGM_n39367_indels_genotypes_effects.csv", header = TRUE)
colnames(effects) = c("VarID", "Effect")

sample_name_gene_name = fread("/Users/sy3115/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/rAF_scripts_and_inputs/PublicationReady_Data/IGM/2023-03-24_11-47-14_IGM_n39367_indels_genotypes_selectcols.csv") 
colnames(sample_name_gene_name) = c("VarID", "geneName", "sampleName", "coveredCtrl", "AC")


annotations = distinct(fread ("/Users/sy3115/Documents/Data/rAF_paper/finalized_project_data/final_dev/2023-04-21_IGM_genename_gnomadpli_gnomadloeuf_omimdisease.csv")) #18231 
colnames(annotations) = c("geneName", "pLI", "oe_lof_upper", "OMIM_disease")

annotations[geneName == "'HTT'"]$OMIM_disease = "Huntington disease, 143100 (3), Autosomal dominant"
annotations[geneName == "'GLTSCR1'"]$geneName = "BICRA"
annotations[geneName == "'FAM46A'"]$geneName = "TENT5A"

merge1_df1 = left_join(input_df1, clinvar, by="VarID")
merge2_df1 = merge(merge1_df1, effects, by = "VarID")
merge3_df1 = merge(merge2_df1, sample_name_gene_name, by = "VarID")
merge4_df1 = merge(merge3_df1, annotations, by= "geneName")

merge1_df2 = left_join(input_df2, clinvar, by="VarID")
merge2_df2 = merge(merge1_df2, effects, by = "VarID")
merge3_df2 = merge(merge2_df2, sample_name_gene_name, by = "VarID")
merge4_df2 = merge(merge3_df2, annotations, by= "geneName")

indels_df1 = merge4_df1
indels_df2 = merge4_df2

## From that list previously, identify a varID of interest.
VarID_of_interest = "4-1980558-GC-G"

## identify the variant of interest within the indels dataframe 
indel_of_interest_df1 = indels_df1[VarID == VarID_of_interest]
indel_of_interest_df2 = indels_df2[VarID == VarID_of_interest]


for (i in bp_range) {
  # get groupID 
  GID_var = paste0 ("GID_bp", i)
  temp_GID = paste0("indel_GID_df_1_bp", i)
    
  assign(temp_GID, indel_of_interest_df1[[GID_var]])
    
  ## Find the indels within same GIDs
  indels_in_GID = paste0 ("indels_in_GID_df_1_bp", i)
    
  assign(indels_in_GID, indels_df1 %>% filter(get(GID_var) == get(temp_GID)))
}

for (i in bp_range) {
  # get groupID 
  GID_var = paste0 ("GID_bp", i)
  temp_GID = paste0("indel_GID_df_2_bp", i)
  
  assign(temp_GID, indel_of_interest_df2[[GID_var]])
  
  ## Find the indels within same GIDs
  indels_in_GID = paste0 ("indels_in_GID_df_2_bp", i)
  
  assign(indels_in_GID, indels_df2 %>% filter(get(GID_var) == get(temp_GID)))
}


view10_df1 = distinct(indels_in_GID_df_1_bp10 %>% select("VarID", "sAF", "rAF_bp10", "V10", "GID_bp10"))
view20_df1 = distinct(indels_in_GID_df_1_bp20 %>% select("VarID", "sAF", "rAF_bp20", "V10", "GID_bp20"))
view30_df1 = distinct(indels_in_GID_df_1_bp30 %>% select("VarID", "sAF", "rAF_bp30", "V10", "GID_bp30"))
view40_df1 = distinct(indels_in_GID_df_1_bp40 %>% select("VarID", "sAF", "rAF_bp40", "V10", "GID_bp40"))

view10_df2 = distinct(indels_in_GID_df_2_bp10 %>% select("VarID", "sAF", "rAF_bp10", "V10", "GID_bp10"))
view20_df2 = distinct(indels_in_GID_df_2_bp20 %>% select("VarID", "sAF", "rAF_bp20", "V10", "GID_bp20"))
view30_df2 = distinct(indels_in_GID_df_2_bp30 %>% select("VarID", "sAF", "rAF_bp30", "V10", "GID_bp30"))
view40_df2 = distinct(indels_in_GID_df_2_bp40 %>% select("VarID", "sAF", "rAF_bp40", "V10", "GID_bp40"))









