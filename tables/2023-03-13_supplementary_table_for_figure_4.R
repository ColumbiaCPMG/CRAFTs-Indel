## Sandy Yang 
## Supplementary Table for Figure 4 
## March 13, 2023 

library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)

gnomad_raw = fread("2023-03-06_figure4_gnomad_input.csv")
IGM_raw = fread("2023-03-06_figure4_IGM_input.csv")

gnomad = gnomad_raw
IGM = IGM_raw

gnomad$rAF5 = round(gnomad_raw$rAF5, 6 )
IGM$rAF5 = round(IGM_raw$rAF5, 6)

## supplementary figure 
## gnomad
gnomad_freq_B_LB = gnomad %>% filter(clinvar_class == "B/LB") %>% group_by(rAF5) %>% summarise(count_B_LB=n())
gnomad_freq_P_LP = gnomad %>% filter(clinvar_class == "P/LP") %>% group_by(rAF5) %>% summarise(count_P_LP=n())
gnomad_freq = merge(gnomad_freq_B_LB, gnomad_freq_P_LP, by = "rAF5", all = TRUE)
# Set NAs to 0 
gnomad_freq[is.na(gnomad_freq)] <- 0 

## IGM
IGM_freq_B_LB = IGM %>% filter(clinvar_class == "B/LB") %>% group_by(rAF5) %>% summarise(count_B_LB=n())
IGM_freq_P_LP = IGM %>% filter(clinvar_class == "P/LP") %>% group_by(rAF5) %>% summarise(count_P_LP=n())
IGM_freq = merge(IGM_freq_B_LB, IGM_freq_P_LP, by = "rAF5", all = TRUE)
# Set NAs to 0 
IGM_freq[is.na(IGM_freq)] <- 0 


## output the supplementary files 
fwrite(gnomad_freq, paste0( Sys.Date(), "_figure4_gnomad_supplementary_table.csv"))
fwrite(IGM_freq, paste0( Sys.Date(), "_figure4_IGM_supplementary_table.csv"))
