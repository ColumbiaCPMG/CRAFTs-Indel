## Sandy Yang 
## May 18, 2023 

###########################################################################################
###########################################################################################
###########################################################################################

## This script is to visualize the rAF_hi indels among chromosomes 

library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
library("ggplot2")
install.packages("ggpubr")
library(ggpubr)

###########################################################################################
###########################################################################################
###########################################################################################

setwd("/Users/User_1/Desktop/rAF_project_2/")

## Read in files 
# rAF_hi_indels_10bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
# rAF_hi_indels_20bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
# rAF_hi_indels_30bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
# rAF_hi_indels_40bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")
# 
# rAF_hi_indels_10bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
# rAF_hi_indels_20bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
# rAF_hi_indels_30bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
# rAF_hi_indels_40bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

rAF_hi_indels_10bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_20bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_30bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_40bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

# df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_lt50bp.csv")
# df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_lt50bp.csv")
df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_lt50bp.csv")

## split varID 
# rAF_hi_indels_10bp_df_1 = rAF_hi_indels_10bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# rAF_hi_indels_20bp_df_1 = rAF_hi_indels_20bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# rAF_hi_indels_30bp_df_1 = rAF_hi_indels_30bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# rAF_hi_indels_40bp_df_1 = rAF_hi_indels_40bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# # 
# rAF_hi_indels_10bp_df_2 = rAF_hi_indels_10bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# rAF_hi_indels_20bp_df_2 = rAF_hi_indels_20bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# rAF_hi_indels_30bp_df_2 = rAF_hi_indels_30bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# rAF_hi_indels_40bp_df_2 = rAF_hi_indels_40bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# # 

rAF_hi_indels_10bp_df_3 = rAF_hi_indels_10bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_20bp_df_3 = rAF_hi_indels_20bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_30bp_df_3 = rAF_hi_indels_30bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_40bp_df_3 = rAF_hi_indels_40bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# df_1 = df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# df_2 = df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
df_3 = df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))

# set bp ranges
bp_range = c("10", "20", "30", "40")

# df_name_1="gnomAD"
# df_name_2="IGM"
df_name_3="UK.BB"

# chromosomes to analyze 
chroms = c("X", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "Y")

## make chroms an ordered factor 
chroms <- factor (chroms, levels = chroms)

rAF_lo = 10^-4 

###########################################################################################
###########################################################################################
###########################################################################################

for (i in bp_range) {
  ## get rAF_hi for that bp range data frame 
  # df_1_rAF_hi = paste0("rAF_hi_indels_", i, "bp_df_1")
  # df_2_rAF_hi = paste0("rAF_hi_indels_", i, "bp_df_2")
  df_3_rAF_hi = paste0("rAF_hi_indels_", i, "bp_df_3")
  
  
  ## filter for rAF_hi indels within each chromosome 
  for (j in chroms) {
    # df_1_rAF_hi_chrom = paste0("df_1_rAF_hi_chrom_", j, "_bp_", i)
    # df_2_rAF_hi_chrom = paste0("df_2_rAF_hi_chrom_", j, "_bp_", i)
    df_3_rAF_hi_chrom = paste0("df_3_rAF_hi_chrom_", j, "_bp_", i)
    
    # assign(df_1_rAF_hi_chrom, nrow(get(df_1_rAF_hi) %>% filter(CHR == j)))
    # assign(df_2_rAF_hi_chrom, nrow(get(df_2_rAF_hi) %>% filter(CHR == j)))
    assign(df_3_rAF_hi_chrom, nrow(get(df_3_rAF_hi) %>% filter(CHR == j)))
    
    ## filter for rAF_lo sAF indels within each chromosome 
    # df_1_rAF_lo_chrom = paste0("df_1_rAF_lo_chrom_", j, "_bp_", i)
    # df_2_rAF_lo_chrom = paste0("df_2_rAF_lo_chrom_", j, "_bp_", i)
    df_3_rAF_lo_chrom = paste0("df_3_rAF_lo_chrom_", j, "_bp_", i)
    
    # assign(df_1_rAF_lo_chrom, nrow(df_1 %>% filter(sAF <= rAF_lo) %>% filter(CHR == j)))
    # assign(df_2_rAF_lo_chrom, nrow(df_2 %>% filter(sAF <= rAF_lo) %>% filter(CHR == j)))
    assign(df_3_rAF_lo_chrom, nrow(df_3 %>% filter(sAF <= rAF_lo) %>% filter(CHR == j)))
    
    ## find percentage of rAF_hi/rAF_lo * 100
    # df_1_prct = paste0("df_1_prct_chrom", j, "_bp_", i)
    # df_2_prct = paste0("df_2_prct_chrom", j, "_bp_", i)
    df_3_prct = paste0("df_3_prct_chrom", j, "_bp_", i)
    
    # assign(df_1_prct, get(df_1_rAF_hi_chrom)/get(df_1_rAF_lo_chrom) * 100 )
    # assign(df_2_prct, get(df_2_rAF_hi_chrom)/get(df_2_rAF_lo_chrom) * 100 )
    assign(df_3_prct, get(df_3_rAF_hi_chrom)/get(df_3_rAF_lo_chrom) * 100 )
    
    
  }
  
  ## make a list for graph data frame 
  # df_1_list = paste0("df_1_list_bp", i)
  # df_2_list = paste0("df_2_list_bp", i)
  df_3_list = paste0("df_3_list_bp", i)
  
  # assign(df_1_list, c())
  # assign(df_2_list, c())
  assign(df_3_list, c())
  
  for (k in chroms) {
    ## get variable 
    # get_temp_df1 = paste0("df_1_prct_chrom", k, "_bp_", i)
    # get_temp_df2 = paste0("df_2_prct_chrom", k, "_bp_", i)
    get_temp_df3 = paste0("df_3_prct_chrom", k, "_bp_", i)
    
    # assign(df_1_list, c(get(df_1_list), get(get_temp_df1)))
    # assign(df_2_list, c(get(df_2_list), get(get_temp_df2)))
    assign(df_3_list, c(get(df_3_list), get(get_temp_df3)))
  }
}

###########################################################################################
###########################################################################################
###########################################################################################

rAF_hi_chr_distr_df1 = data.frame(x = chroms, y1 = df_1_list_bp10, y2 = df_1_list_bp20, y3 = df_1_list_bp30, y4 = df_1_list_bp40)

rAF_hi_chr_distr_df2 = data.frame(x = chroms, y1 = df_2_list_bp10, y2 = df_2_list_bp20, y3 = df_2_list_bp30, y4 = df_2_list_bp40)

rAF_hi_chr_distr_df3 = data.frame(x = chroms, y1 = df_3_list_bp10, y2 = df_3_list_bp20, y3 = df_3_list_bp30, y4 = df_3_list_bp40)


rAF_hi_chr_distr_plot1 = ggplot(rAF_hi_chr_distr_df1, aes(x)) +
  geom_line(aes(y = y1, color = "10bp region"), group = 1) +
  geom_line(aes(y = y2, color = "20bp region"), group = 1) +
  geom_line(aes(y = y3, color = "30bp region"), group = 1) +
  geom_line(aes(y = y4, color = "40bp region"), group = 1) +
  ylim (0, 30) +
  labs (title = str_wrap(paste0("Percentage of rAF_lo Indels that are rAF_hi Across Chromosomes in the ", df_name_1, " dataset"), width = 60),
        x = "Chromosomes",
        y = "Percentage")

rAF_hi_chr_distr_plot2 = ggplot(rAF_hi_chr_distr_df2, aes(x)) +
  geom_line(aes(y = y1, color = "10bp region"), group = 1) +
  geom_line(aes(y = y2, color = "20bp region"), group = 1) +
  geom_line(aes(y = y3, color = "30bp region"), group = 1) +
  geom_line(aes(y = y4, color = "40bp region"), group = 1) +
  ylim (0,30) +
  labs (title = str_wrap(paste0("Percentage of rAF_lo Indels that are rAF_hi Across Chromosomes in the ", df_name_2, " dataset"), width = 60),
        x = "Chromosomes",
        y = "Percentage")

rAF_hi_chr_distr_plot3 = ggplot(rAF_hi_chr_distr_df3, aes(x)) +
  geom_line(aes(y = y1, color = "10bp region"), group = 1) +
  geom_line(aes(y = y2, color = "20bp region"), group = 1) +
  geom_line(aes(y = y3, color = "30bp region"), group = 1) +
  geom_line(aes(y = y4, color = "40bp region"), group = 1) +
  ylim (0,30) +
  labs (title = str_wrap(paste0("Percentage of rAF_lo Indels that are rAF_hi Across Chromosomes in the ", df_name_3, " dataset"), width = 60),
        x = "Chromosomes",
        y = "Percentage")





###########################################################################################
###########################################################################################
###########################################################################################
#Create the Data for the merged plot 

rAF_hi_chr_distr_df1$Plot <- "gnomAD"
rAF_hi_chr_distr_df2$Plot <- "IGM"
rAF_hi_chr_distr_df3$Plot <- "UK.BB"


combined_data <- rbind(rAF_hi_chr_distr_df1, rAF_hi_chr_distr_df2, rAF_hi_chr_distr_df3)


combined_data_long <- combined_data %>%
  gather(key = "variable", value = "value", -x, -Plot)



combined_plot <- ggplot(combined_data_long, aes(x, value, color = variable, group = variable)) +
  geom_line() +
  facet_grid(Plot ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Chromosome", y = "Percentage", title = "Percentage of rAF_lo Indels that are rAF_hi Across Chromosomes in the gnomAD, IGM, and UK.BB Data Sets.") +
  theme_minimal()


ggsave("FigS1.jpg", width = 80, height = 60, units = c("cm"), dpi = 300)
