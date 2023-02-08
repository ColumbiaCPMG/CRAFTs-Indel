## Author: Sandy Yang
## December 6, 2022
## This code is the input for figure 2 of the indel paper 
## Here, we want to find the number of indels with the gnomAD and IGM data sets that are rare by only standard AF (sAF) or those that are rare by both sAF and rAF (regional allele frequency)

library(dplyr)

### define rare
### rare is 1 * 10^-4 for gnomAD too because we multiplied the sAF/rAF of the dataset by 2 to account for the AN value discrepancy 
rare = (1 * 10^-4)

gnomad_indels = read.csv("2022-12-30_gnomAD_sAF_rAF.csv", sep=',', header = TRUE)
IGM_indels = read.csv("2022-12-30_IGM_sAF_rAF.csv", sep=',', header = TRUE)

## total number of indels rare by sAF 
tot_rare_sAF_gnomad = nrow(gnomad_indels %>% filter(sAF <= rare))
tot_rare_sAF_IGM = nrow(IGM_indels %>% filter( sAF <= rare))

## Find the number of indels rare by sAF and common by rAF 
## gnomAD
rare_sAF_common_rAF_bp5 = nrow(gnomad_indels %>% filter(sAF <= rare & rAF5 > rare))
rare_sAF_common_rAF_bp10 = nrow(gnomad_indels %>% filter(sAF <= rare & rAF10 > rare))
rare_sAF_common_rAF_bp15 = nrow(gnomad_indels %>% filter(sAF <= rare & rAF15 > rare))
rare_sAF_common_rAF_bp20 = nrow(gnomad_indels %>% filter(sAF <= rare & rAF20 > rare))

percent_rare_sAF_common_rAF_bp5 = rare_sAF_common_rAF_bp5/tot_rare_sAF_gnomad * 100
percent_rare_sAF_common_rAF_bp10 = rare_sAF_common_rAF_bp10/tot_rare_sAF_gnomad * 100
percent_rare_sAF_common_rAF_bp15 = rare_sAF_common_rAF_bp15/tot_rare_sAF_gnomad * 100
percent_rare_sAF_common_rAF_bp20 = rare_sAF_common_rAF_bp20/tot_rare_sAF_gnomad * 100

## IGM
rare_sAF_common_rAF_bp5_IGM = nrow(IGM_indels %>% filter(sAF <= rare & rAF5 > rare))
rare_sAF_common_rAF_bp10_IGM = nrow(IGM_indels %>% filter(sAF <= rare & rAF10 > rare))
rare_sAF_common_rAF_bp15_IGM = nrow(IGM_indels %>% filter(sAF <= rare & rAF15 > rare))
rare_sAF_common_rAF_bp20_IGM = nrow(IGM_indels %>% filter(sAF <= rare & rAF20 > rare))

percent_rare_sAF_common_rAF_bp5_IGM = rare_sAF_common_rAF_bp5_IGM/tot_rare_sAF_IGM * 100
percent_rare_sAF_common_rAF_bp10_IGM = rare_sAF_common_rAF_bp10_IGM/tot_rare_sAF_IGM * 100
percent_rare_sAF_common_rAF_bp15_IGM = rare_sAF_common_rAF_bp15_IGM/tot_rare_sAF_IGM * 100
percent_rare_sAF_common_rAF_bp20_IGM = rare_sAF_common_rAF_bp20_IGM/tot_rare_sAF_IGM * 100

## Find the number of indels rare by both sAF and rAF 
## gnomAD
rare_sAF_rare_rAF_bp5 = nrow(gnomad_indels %>% filter(sAF <= rare & rAF5 <= rare ))
rare_sAF_rare_rAF_bp10 = nrow(gnomad_indels %>% filter(sAF <= rare & rAF10 <= rare ))
rare_sAF_rare_rAF_bp15 = nrow(gnomad_indels %>% filter(sAF <= rare & rAF15 <= rare ))
rare_sAF_rare_rAF_bp20 = nrow(gnomad_indels %>% filter(sAF <= rare & rAF20 <= rare ))

percent_rare_sAF_rare_rAF_bp5 = rare_sAF_rare_rAF_bp5 / tot_rare_sAF_gnomad * 100
percent_rare_sAF_rare_rAF_bp10 = rare_sAF_rare_rAF_bp10 / tot_rare_sAF_gnomad * 100
percent_rare_sAF_rare_rAF_bp15 = rare_sAF_rare_rAF_bp15 / tot_rare_sAF_gnomad * 100
percent_rare_sAF_rare_rAF_bp20 = rare_sAF_rare_rAF_bp20 / tot_rare_sAF_gnomad * 100

## IGM 
rare_sAF_rare_rAF_bp5_IGM = nrow(IGM_indels %>% filter(sAF <= rare & rAF5 <= rare ))
rare_sAF_rare_rAF_bp10_IGM = nrow(IGM_indels %>% filter(sAF <= rare & rAF10 <= rare ))
rare_sAF_rare_rAF_bp15_IGM = nrow(IGM_indels %>% filter(sAF <= rare & rAF15 <= rare ))
rare_sAF_rare_rAF_bp20_IGM = nrow(IGM_indels %>% filter(sAF <= rare & rAF20 <= rare ))

percent_rare_sAF_rare_rAF_bp5_IGM = rare_sAF_rare_rAF_bp5_IGM / tot_rare_sAF_IGM * 100
percent_rare_sAF_rare_rAF_bp10_IGM = rare_sAF_rare_rAF_bp10_IGM / tot_rare_sAF_IGM * 100
percent_rare_sAF_rare_rAF_bp15_IGM = rare_sAF_rare_rAF_bp15_IGM / tot_rare_sAF_IGM * 100
percent_rare_sAF_rare_rAF_bp20_IGM = rare_sAF_rare_rAF_bp20_IGM / tot_rare_sAF_IGM * 100

## sanity checks: should add up to the total number of rare sAF 
## GNOMAD
percent_rare_sAF_common_rAF_bp5 + percent_rare_sAF_rare_rAF_bp5
percent_rare_sAF_common_rAF_bp10 + percent_rare_sAF_rare_rAF_bp10
percent_rare_sAF_common_rAF_bp15 + percent_rare_sAF_rare_rAF_bp15
percent_rare_sAF_common_rAF_bp20 + percent_rare_sAF_rare_rAF_bp20

tot_rare_sAF_gnomad
rare_sAF_common_rAF_bp5 + rare_sAF_rare_rAF_bp5
rare_sAF_common_rAF_bp10 + rare_sAF_rare_rAF_bp10
rare_sAF_common_rAF_bp15 + rare_sAF_rare_rAF_bp15
rare_sAF_common_rAF_bp20 + rare_sAF_rare_rAF_bp20

## IGM 
percent_rare_sAF_rare_rAF_bp5_IGM + percent_rare_sAF_common_rAF_bp5_IGM
percent_rare_sAF_rare_rAF_bp10_IGM + percent_rare_sAF_common_rAF_bp10_IGM
percent_rare_sAF_rare_rAF_bp15_IGM + percent_rare_sAF_common_rAF_bp15_IGM
percent_rare_sAF_rare_rAF_bp20_IGM + percent_rare_sAF_common_rAF_bp20_IGM

tot_rare_sAF_IGM
rare_sAF_rare_rAF_bp5_IGM + rare_sAF_common_rAF_bp5_IGM
rare_sAF_rare_rAF_bp10_IGM + rare_sAF_common_rAF_bp10_IGM
rare_sAF_rare_rAF_bp15_IGM + rare_sAF_common_rAF_bp15_IGM
rare_sAF_rare_rAF_bp20_IGM + rare_sAF_common_rAF_bp20_IGM

## output csv files 
## gnomAD 
windows = c("window10", "window20", "window30", "window40")

percent_rare_sAF_rare_rAF_gnomad = c(percent_rare_sAF_rare_rAF_bp5, percent_rare_sAF_rare_rAF_bp10, percent_rare_sAF_rare_rAF_bp15, percent_rare_sAF_rare_rAF_bp20)

percent_rare_sAF_common_rAF_gnomad = c(percent_rare_sAF_common_rAF_bp5, percent_rare_sAF_common_rAF_bp10, percent_rare_sAF_common_rAF_bp15, percent_rare_sAF_common_rAF_bp20)

rare_sAF_rare_rAF_gnomad = c(rare_sAF_rare_rAF_bp5, rare_sAF_rare_rAF_bp10, rare_sAF_rare_rAF_bp15, rare_sAF_rare_rAF_bp20)

rare_sAF_common_rAF_gnomad = c(rare_sAF_common_rAF_bp5, rare_sAF_common_rAF_bp10, rare_sAF_common_rAF_bp15, rare_sAF_common_rAF_bp20)

## dataframe for gnomad 
gnomad_output_fig_2 = data.frame(windows, percent_rare_sAF_rare_rAF_gnomad, percent_rare_sAF_common_rAF_gnomad, rare_sAF_rare_rAF_gnomad, rare_sAF_common_rAF_gnomad)

write.table(gnomad_output_fig_2, paste0( Sys.Date(), "_gnomad_fig2_input.csv"), row.names = FALSE, sep = "," )

percent_rare_sAF_rare_rAF_IGM = c(percent_rare_sAF_rare_rAF_bp5_IGM, percent_rare_sAF_rare_rAF_bp10_IGM, percent_rare_sAF_rare_rAF_bp15_IGM, percent_rare_sAF_rare_rAF_bp20_IGM)

percent_rare_sAF_common_rAF_IGM = c(percent_rare_sAF_common_rAF_bp5_IGM, percent_rare_sAF_common_rAF_bp10_IGM, percent_rare_sAF_common_rAF_bp15_IGM, percent_rare_sAF_common_rAF_bp20_IGM)

rare_sAF_rare_rAF_IGM = c(rare_sAF_rare_rAF_bp5_IGM, rare_sAF_rare_rAF_bp10_IGM, rare_sAF_rare_rAF_bp15_IGM, rare_sAF_rare_rAF_bp20_IGM)

rare_sAF_common_rAF_IGM = c(rare_sAF_common_rAF_bp5_IGM, rare_sAF_common_rAF_bp10_IGM, rare_sAF_common_rAF_bp15_IGM, rare_sAF_common_rAF_bp20_IGM)

## dataframe for IGM
IGM_output_fig_2 = data.frame(windows, percent_rare_sAF_rare_rAF_IGM, percent_rare_sAF_common_rAF_IGM, rare_sAF_rare_rAF_IGM, rare_sAF_common_rAF_IGM)

write.table(IGM_output_fig_2, paste0( Sys.Date(), "_IGM_fig2_input.csv"), row.names = FALSE, sep = "," )







  
  