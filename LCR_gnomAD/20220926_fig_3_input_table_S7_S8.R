## Author: Sandy Yang
## January 6, 2023 
## This code is to determine the percentage of gnomad indels that are rare by sAC/sAF and common by rAC/rAF, that are within the LCR region 
remove.packages("rlang")
install.packages("rlang")
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)

## these are manual inputs, but these numbers can also be found by finding the rows of these files below
#rare_sAF_common_rAF_LCR_bp5 = 28314 
#rare_sAF_common_rAF_LCR_bp10 = 32588
#rare_sAF_common_rAF_LCR_bp15 = 34518
#rare_sAF_common_rAF_LCR_bp20 = 35112

rare_sAF_common_rAF_LCR_bp5 = nrow(fread("2023-03-06_LCR_rare_sAF_common_rAF_5.bed", sep = ",", header = FALSE))
rare_sAF_common_rAF_LCR_bp10 = nrow(fread("2023-03-06_LCR_rare_sAF_common_rAF_10.bed", sep = ",", header = FALSE))
rare_sAF_common_rAF_LCR_bp15 = nrow(fread("2023-03-06_LCR_rare_sAF_common_rAF_15.bed", sep = ",", header = FALSE))
rare_sAF_common_rAF_LCR_bp20 = nrow(fread("2023-03-06_LCR_rare_sAF_common_rAF_20.bed", sep = ",", header = FALSE))

## Since the rAF and sAF are both multiplied by two, the rare threshold is also multiplied by 2 
rare <- (1 * 10^-4)

#gnomad_indels = read.csv("/Users/sy3115/Documents/Data/20220926_gnomad_LCR_bedtools_intersect_official/20220926_gnomad_srAF_double_rare_threshold_10_4.csv", sep=',', header = TRUE)

gnomad_indels = read.csv("2023-03-06_gnomad_indels_AF.csv", sep=',', header = TRUE)

total_rare_sAF_common_rAF_5 = nrow(gnomad_indels %>% filter (sAF <= rare & rAF5 > rare))
total_rare_sAF_common_rAF_10 = nrow(gnomad_indels %>% filter (sAF <= rare & rAF10 > rare))
total_rare_sAF_common_rAF_15 = nrow(gnomad_indels %>% filter (sAF <= rare & rAF15 > rare))
total_rare_sAF_common_rAF_20 = nrow(gnomad_indels %>% filter (sAF <= rare & rAF20 > rare))
  
percent_rare_sAF_common_rAF_LCR_bp5 = (rare_sAF_common_rAF_LCR_bp5/total_rare_sAF_common_rAF_5) *100

percent_rare_sAF_common_rAF_LCR_bp10 = (rare_sAF_common_rAF_LCR_bp10/total_rare_sAF_common_rAF_10)*100
  
percent_rare_sAF_common_rAF_LCR_bp15 = (rare_sAF_common_rAF_LCR_bp15/total_rare_sAF_common_rAF_15)*100
  
percent_rare_sAF_common_rAF_LCR_bp20 = (rare_sAF_common_rAF_LCR_bp20/total_rare_sAF_common_rAF_20)*100

percent_rare_sAF_common_rAF_LCR_bp5
percent_rare_sAF_common_rAF_LCR_bp10
percent_rare_sAF_common_rAF_LCR_bp15
percent_rare_sAF_common_rAF_LCR_bp20


## Find rare sAF and common rAF indels not in LCR
not_in_LCR_rare_sAF_common_rAF_5 = total_rare_sAF_common_rAF_5 - rare_sAF_common_rAF_LCR_bp5
not_in_LCR_rare_sAF_common_rAF_10 = total_rare_sAF_common_rAF_10 - rare_sAF_common_rAF_LCR_bp10
not_in_LCR_rare_sAF_common_rAF_15 = total_rare_sAF_common_rAF_15 - rare_sAF_common_rAF_LCR_bp15
not_in_LCR_rare_sAF_common_rAF_20 = total_rare_sAF_common_rAF_20 - rare_sAF_common_rAF_LCR_bp20

not_in_LCR_rare_sAF_common_rAF_5
not_in_LCR_rare_sAF_common_rAF_10
not_in_LCR_rare_sAF_common_rAF_15
not_in_LCR_rare_sAF_common_rAF_20

percent_not_in_LCR_5 = (not_in_LCR_rare_sAF_common_rAF_5 / total_rare_sAF_common_rAF_5) * 100
percent_not_in_LCR_10 = (not_in_LCR_rare_sAF_common_rAF_10 / total_rare_sAF_common_rAF_10) * 100
percent_not_in_LCR_15 = (not_in_LCR_rare_sAF_common_rAF_15 / total_rare_sAF_common_rAF_15) * 100
percent_not_in_LCR_20 = (not_in_LCR_rare_sAF_common_rAF_20 / total_rare_sAF_common_rAF_20) * 100

percent_not_in_LCR_5
percent_not_in_LCR_10 
percent_not_in_LCR_15
percent_not_in_LCR_20

# sanity check: should add up to 100%
percent_rare_sAF_common_rAF_LCR_bp5 + percent_not_in_LCR_5
percent_rare_sAF_common_rAF_LCR_bp10 + percent_not_in_LCR_10
percent_rare_sAF_common_rAF_LCR_bp15 + percent_not_in_LCR_15
percent_rare_sAF_common_rAF_LCR_bp20 + percent_not_in_LCR_20

not_in_LCR_rare_sAF_common_rAF_5 + rare_sAF_common_rAF_LCR_bp5
total_rare_sAF_common_rAF_5

not_in_LCR_rare_sAF_common_rAF_10 + rare_sAF_common_rAF_LCR_bp10
total_rare_sAF_common_rAF_10

not_in_LCR_rare_sAF_common_rAF_15 + rare_sAF_common_rAF_LCR_bp15
total_rare_sAF_common_rAF_15

not_in_LCR_rare_sAF_common_rAF_20 + rare_sAF_common_rAF_LCR_bp20
total_rare_sAF_common_rAF_20

### output a csv file as an input for the side-by-side graph for LCR (figure 4)
window = c("window10", "window20", "window30", "window40")

percent_in_LCR = c(percent_rare_sAF_common_rAF_LCR_bp5, percent_rare_sAF_common_rAF_LCR_bp10, percent_rare_sAF_common_rAF_LCR_bp15, percent_rare_sAF_common_rAF_LCR_bp20)

percent_outside_LCR = c(percent_not_in_LCR_5, percent_not_in_LCR_10, percent_not_in_LCR_15, percent_not_in_LCR_20)

in_LCR = c(rare_sAF_common_rAF_LCR_bp5, rare_sAF_common_rAF_LCR_bp10, rare_sAF_common_rAF_LCR_bp15, rare_sAF_common_rAF_LCR_bp20)

outside_LCR = c(not_in_LCR_rare_sAF_common_rAF_5, not_in_LCR_rare_sAF_common_rAF_10,not_in_LCR_rare_sAF_common_rAF_15, not_in_LCR_rare_sAF_common_rAF_20)

gnomad_output_df = data.frame(window, percent_in_LCR, percent_outside_LCR, in_LCR, outside_LCR)

write.table(gnomad_output_df, paste0( Sys.Date(), "_gnomad_fig3_input.csv"), row.names = FALSE, sep = "," )

  
######## More tests ########
# total_rare_sAF_common_rAF_5 / nrow(gnomad_indels) * 100 
# total_rare_sAF_common_rAF_10 / nrow(gnomad_indels) * 100
# total_rare_sAF_common_rAF_15 / nrow(gnomad_indels) * 100
# total_rare_sAF_common_rAF_20 / nrow(gnomad_indels) * 100
# 
# ## find indels rare by both sAF and rAF 
# total_rare_sAF_rare_rAF_5 = nrow(gnomad_indels %>% filter (sAF <= rare & rAF5 <= rare))
# total_rare_sAF_rare_rAF_10 = nrow(gnomad_indels %>% filter (sAF <= rare & rAF10 <= rare))
# total_rare_sAF_rare_rAF_15 = nrow(gnomad_indels %>% filter (sAF <= rare & rAF15 <= rare))
# total_rare_sAF_rare_rAF_20 = nrow(gnomad_indels %>% filter (sAF <= rare & rAF20 <= rare))
# 
# ## find indels that are rare by sAF 
# rare_sAF =  nrow(gnomad_indels %>% filter (sAF <= rare))
# 
# ## see if indels are rare by rAF5 but common by sAF
# ## yes! This is possible because AF calculation uses a mean AN rather than absolute AN for position 
# rare_rAF5_common_sAF = gnomad_indels %>% filter(sAF > rare & rAF5 <= rare)

