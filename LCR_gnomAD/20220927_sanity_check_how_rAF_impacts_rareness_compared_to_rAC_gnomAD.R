## Author: Sandy Yang 
## September 27, 2022 
## This code is to check the percentage of indels that are rare based on rAC and common based on rAF
getwd();
setwd("~/Desktop/");
getwd();
## header = true means that the first row of the table includes the header and therefore it will not be counted when we do nrow
gnomad_indels_with_rAF_formatted = read.csv(file = '20220926_gnomad_srAF_double_rare_threshold_10_4.csv', sep = ',', header = TRUE)
nrow (gnomad_indels_with_rAF_formatted)
## fix the csv file if necessary 
##gnomad_indels_with_rAF_formatted = subset (gnomad_indels_with_rAF, select = -c(VarID, X))

## install correct packages for filter
remove.packages("rlang")
install.packages("rlang")
library(rlang)
library(dplyr)

rare_AC = 12 
rare_AF = 1 * 10^-4

## Filter data for those that are rare by rAC5
gnomad_indels_rare_by_rAC_5 <- gnomad_indels_with_rAF_formatted %>% filter(AlleleCount_bp5<=rare_AC)
percentage_rare_by_rAC_5 = nrow(gnomad_indels_rare_by_rAC_5) / nrow(gnomad_indels_with_rAF_formatted) * 100
percentage_rare_by_rAC_5

## Filter data for those that are rare by rAF5
gnomad_indels_rare_by_rAF5 <- gnomad_indels_with_rAF_formatted %>% filter(rAF5<= rare_AF)
percentage_rare_by_rAF5 = nrow(gnomad_indels_rare_by_rAF5) / nrow(gnomad_indels_with_rAF_formatted) * 100
percentage_rare_by_rAF5

## Filter data for those that are rare by rAC10
gnomad_indels_rare_by_rAC_10 <- gnomad_indels_with_rAF_formatted %>% filter(AlleleCount_bp10<=rare_AC)
percentage_rare_by_rAC_10 = nrow(gnomad_indels_rare_by_rAC_10) / nrow(gnomad_indels_with_rAF_formatted) * 100
percentage_rare_by_rAC_10

## Filter data for those that are rare by rAF10
gnomad_indels_rare_by_rAF10 <- gnomad_indels_with_rAF_formatted %>% filter(rAF10<= rare_AF)
percentage_rare_by_rAF10 = nrow(gnomad_indels_rare_by_rAF10) / nrow(gnomad_indels_with_rAF_formatted) * 100
percentage_rare_by_rAF10

## Filter data for those that are rare by rAC15
gnomad_indels_rare_by_rAC_15 <- gnomad_indels_with_rAF_formatted %>% filter(AlleleCount_bp15<=rare_AC)
percentage_rare_by_rAC_15 = nrow(gnomad_indels_rare_by_rAC_15) / nrow(gnomad_indels_with_rAF_formatted) * 100
percentage_rare_by_rAC_15

## Filter data for those that are rare by rAF15
gnomad_indels_rare_by_rAF15 <- gnomad_indels_with_rAF_formatted %>% filter(rAF15<= rare_AF)
percentage_rare_by_rAF15 = nrow(gnomad_indels_rare_by_rAF15) / nrow(gnomad_indels_with_rAF_formatted) * 100
percentage_rare_by_rAF15

## Filter data for those that are rare by rAC20
gnomad_indels_rare_by_rAC_20 <- gnomad_indels_with_rAF_formatted %>% filter(AlleleCount_bp20<=rare_AC)
percentage_rare_by_rAC_20 = nrow(gnomad_indels_rare_by_rAC_20) / nrow(gnomad_indels_with_rAF_formatted) * 100
percentage_rare_by_rAC_20

## Filter data for those that are rare by rAF20
gnomad_indels_rare_by_rAF20 <- gnomad_indels_with_rAF_formatted %>% filter(rAF20<= rare_AF)
percentage_rare_by_rAF20 = nrow(gnomad_indels_rare_by_rAF20) / nrow(gnomad_indels_with_rAF_formatted) * 100
percentage_rare_by_rAF20

## Find the discrepancy between the two dataframes 
## when threshold for sAF is set to 5 * 10^-5, the discrepancy all becomes that indels are being classified as rare sAC but are still common sAF 
## Since the sAF and rAF are doubled, we also double the threshold to 1 * 10^-4
## The sAF becomes more stringent for rareness filter 
## compare discrepancy for bp 5 
rare_rAC5_common_rAF5 <- gnomad_indels_with_rAF_formatted %>% filter (AlleleCount_bp5 <= rare_AC & rAF5 > rare_AF)
nrow(rare_rAC5_common_rAF5)

common_rAC5_rare_sAF5 <- gnomad_indels_with_rAF_formatted %>% filter (AlleleCount_bp5 > rare_AC & rAF5 <= rare_AF)
nrow (common_rAC5_rare_sAF5)

## compare discrepancy for bp 10
rare_rAC10_common_rAF10 <- gnomad_indels_with_rAF_formatted %>% filter (AlleleCount_bp10 <= rare_AC & rAF10 > rare_AF)
nrow(rare_rAC10_common_rAF10)

common_rAC5_rare_rAF10 <- gnomad_indels_with_rAF_formatted %>% filter (AlleleCount_bp10 > rare_AC & rAF10 <= rare_AF)
nrow (common_rAC5_rare_rAF10)

## compare discrepancy for bp 15
rare_rAC15_common_rAF15 <- gnomad_indels_with_rAF_formatted %>% filter (AlleleCount_bp15 <= rare_AC & rAF15 > rare_AF)
nrow(rare_rAC15_common_rAF15)

common_rAC15_rare_sAF15 <- gnomad_indels_with_rAF_formatted %>% filter (AlleleCount_bp15 > rare_AC & rAF15 <= rare_AF)
nrow (common_rAC15_rare_sAF15)

## compare discrepancy for bp 20 
rare_rAC20_common_rAF20 <- gnomad_indels_with_rAF_formatted %>% filter (AlleleCount_bp20 <= rare_AC & rAF20 > rare_AF)
nrow(rare_rAC20_common_rAF20)

common_rAC20_rare_sAF20 <- gnomad_indels_with_rAF_formatted %>% filter (AlleleCount_bp20 > rare_AC & rAF20 <= rare_AF)
nrow (common_rAC20_rare_sAF20)

## Output the dataframes as a csv
write.csv(rare_rAC5_common_rAF5, 'rare_rAC5_common_rAF5_threshold_5_10_5_no_nonzero_ACAN.csv', row.names = FALSE, quote = FALSE)
write.csv(rare_rAC10_common_rAF10, 'rare_rAC10_common_rAF10_threshold_5_10_5_no_nonzero_ACAN.csv', row.names = FALSE, quote = FALSE)
write.csv(rare_rAC15_common_rAF15, 'rare_rAC15_common_rAF15_threshold_5_10_5_no_nonzero_ACAN.csv', row.names = FALSE, quote = FALSE)
write.csv(rare_rAC20_common_rAF20, 'rare_rAC20_common_rAF20_threshold_5_10_5_no_nonzero_ACAN.csv', row.names = FALSE, quote = FALSE)

## Check to see if all the indels in gnomad_indels_rare_by_sAF are in gnomad_indels_rare_by_AC
## concatenate the varID 
gnomad_indels_rare_by_rAC_5$VarID <- paste(gnomad_indels_rare_by_rAC_5$CHR, gnomad_indels_rare_by_rAC_5$POS, gnomad_indels_rare_by_rAC_5$REF, gnomad_indels_rare_by_rAC_5$ALT)
gnomad_indels_rare_by_rAF5$VarID <- paste(gnomad_indels_rare_by_rAF5$CHR, gnomad_indels_rare_by_rAF5$POS, gnomad_indels_rare_by_rAF5$REF, gnomad_indels_rare_by_rAF5$ALT)
## Check to see if the indels in gnomad_indels_rare_by_sAF are in gnomad_indels_rare_by_AC using VarID
indels_in_rare_by_rAF5_but_not_rare_by_rAF5 = subset(gnomad_indels_rare_by_rAF5, !(gnomad_indels_rare_by_rAF5$VarID %in% gnomad_indels_rare_by_rAC_5$VarID))
nrow(indels_in_rare_by_rAF5_but_not_rare_by_rAF5)

## concatenate the varID 
gnomad_indels_rare_by_rAC_10$VarID <- paste(gnomad_indels_rare_by_rAC_10$CHR, gnomad_indels_rare_by_rAC_10$POS, gnomad_indels_rare_by_rAC_10$REF, gnomad_indels_rare_by_rAC_10$ALT)
gnomad_indels_rare_by_rAF10$VarID <- paste(gnomad_indels_rare_by_rAF10$CHR, gnomad_indels_rare_by_rAF10$POS, gnomad_indels_rare_by_rAF10$REF, gnomad_indels_rare_by_rAF10$ALT)
## Check to see if the indels in gnomad_indels_rare_by_sAF are in gnomad_indels_rare_by_AC using VarID
indels_in_rare_by_rAF10_but_not_rare_by_rAF10 = subset(gnomad_indels_rare_by_rAF10, !(gnomad_indels_rare_by_rAF10$VarID %in% gnomad_indels_rare_by_rAC_10$VarID))
nrow(indels_in_rare_by_rAF10_but_not_rare_by_rAF10)

## concatenate the varID 
gnomad_indels_rare_by_rAC_15$VarID <- paste(gnomad_indels_rare_by_rAC_15$CHR, gnomad_indels_rare_by_rAC_15$POS, gnomad_indels_rare_by_rAC_15$REF, gnomad_indels_rare_by_rAC_15$ALT)
gnomad_indels_rare_by_rAF15$VarID <- paste(gnomad_indels_rare_by_rAF15$CHR, gnomad_indels_rare_by_rAF15$POS, gnomad_indels_rare_by_rAF15$REF, gnomad_indels_rare_by_rAF15$ALT)
## Check to see if the indels in gnomad_indels_rare_by_sAF are in gnomad_indels_rare_by_AC using VarID
indels_in_rare_by_rAF15_but_not_rare_by_rAF15 = subset(gnomad_indels_rare_by_rAF15, !(gnomad_indels_rare_by_rAF15$VarID %in% gnomad_indels_rare_by_rAC_15$VarID))
nrow(indels_in_rare_by_rAF15_but_not_rare_by_rAF15)

## concatenate the varID 
gnomad_indels_rare_by_rAC_20$VarID <- paste(gnomad_indels_rare_by_rAC_20$CHR, gnomad_indels_rare_by_rAC_20$POS, gnomad_indels_rare_by_rAC_20$REF, gnomad_indels_rare_by_rAC_20$ALT)
gnomad_indels_rare_by_rAF20$VarID <- paste(gnomad_indels_rare_by_rAF20$CHR, gnomad_indels_rare_by_rAF20$POS, gnomad_indels_rare_by_rAF20$REF, gnomad_indels_rare_by_rAF20$ALT)
## Check to see if the indels in gnomad_indels_rare_by_sAF are in gnomad_indels_rare_by_AC using VarID
indels_in_rare_by_rAF20_but_not_rare_by_rAF20 = subset(gnomad_indels_rare_by_rAF20, !(gnomad_indels_rare_by_rAF20$VarID %in% gnomad_indels_rare_by_rAC_20$VarID))
nrow(indels_in_rare_by_rAF20_but_not_rare_by_rAF20)




