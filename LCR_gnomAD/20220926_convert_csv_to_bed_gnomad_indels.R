## January 6, 2023 
## Author: Sandy Yang
## This code converts gnomad indel csv files into bed files genome-wide

remove.packages("rlang")
install.packages("rlang")
library(rlang)
library(dplyr)


## Since the rAF and sAF are both multiplied by two, the rare threshold is also multiplied by 2 
rare <- (1 * 10^-4)

indels_AF = read.csv(file = '2023-03-06_gnomad_indels_AF.csv', sep = ',', header = TRUE)

## find those that are rare by sAF but common by rAF for bp5
rare_sAF_common_rAF_5 = indels_AF %>% filter(sAF <= rare & rAF5 > rare)
rare_sAF_common_rAF_5_bedfile = subset (rare_sAF_common_rAF_5, select = c(CHR, POS))
rare_sAF_common_rAF_5_bedfile$end_pos = rare_sAF_common_rAF_5_bedfile$POS

write.table(rare_sAF_common_rAF_5_bedfile, "2023-03-06_rare_sAF_common_rAF_5_bedfile.bed", quote = F, row.names = F, col.names = F, sep = "\t")

## find those that are rare by sAF but common by rAF for bp10
rare_sAF_common_rAF_10 = indels_AF %>% filter(sAF <= rare & rAF10 > rare)
rare_sAF_common_rAF_10_bedfile = subset (rare_sAF_common_rAF_10, select = c(CHR, POS))
rare_sAF_common_rAF_10_bedfile$end_pos = rare_sAF_common_rAF_10_bedfile$POS

write.table(rare_sAF_common_rAF_10_bedfile, "2023-03-06_rare_sAF_common_rAF_10_bedfile.bed", quote = F, row.names = F, col.names = F, sep = "\t")

## find those that are rare by sAF but common by rAF for bp15
rare_sAF_common_rAF_15 = indels_AF %>% filter(sAF <= rare & rAF15 > rare)
rare_sAF_common_rAF_15_bedfile = subset (rare_sAF_common_rAF_15, select = c(CHR, POS))
rare_sAF_common_rAF_15_bedfile$end_pos = rare_sAF_common_rAF_15_bedfile$POS

write.table(rare_sAF_common_rAF_15_bedfile, "2023-03-06_rare_sAF_common_rAF_15_bedfile.bed", quote = F, row.names = F, col.names = F, sep = "\t")

## find those that are rare by sAF but common by rAF for bp20 
rare_sAF_common_rAF_20 = indels_AF %>% filter(sAF <= rare & rAF20 > rare)
rare_sAF_common_rAF_20_bedfile = subset (rare_sAF_common_rAF_20, select = c(CHR, POS))
rare_sAF_common_rAF_20_bedfile$end_pos = rare_sAF_common_rAF_20_bedfile$POS

write.table(rare_sAF_common_rAF_20_bedfile, "2023-03-06_rare_sAF_common_rAF_20_bedfile.bed", quote = F, row.names = F, col.names = F, sep = "\t")


