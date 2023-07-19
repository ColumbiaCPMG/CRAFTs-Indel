## Sandy Yang 
## May 18, 2023 

###########################################################################################
###########################################################################################
###########################################################################################

## This script is to visualize the suspicious indels among chromosomes 

library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library("ggplot2")
install.packages("ggpubr")
library(ggpubr)

###########################################################################################
###########################################################################################
###########################################################################################

setwd("")

## Read in files 
#suspicious_indels_10bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_SuspiciousIndels.lt50bp.csv")
#suspicious_indels_20bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_SuspiciousIndels.lt50bp.csv")
#suspicious_indels_30bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_SuspiciousIndels.lt50bp.csv")
#suspicious_indels_40bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_SuspiciousIndels.lt50bp.csv")

#suspicious_indels_10bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_SuspiciousIndels.lt50bp.csv")
#suspicious_indels_20bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp20_SuspiciousIndels.lt50bp.csv")
#suspicious_indels_30bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp30_SuspiciousIndels.lt50bp.csv")
#suspicious_indels_40bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp40_SuspiciousIndels.lt50bp.csv")

suspicious_indels_10bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_SuspiciousIndels.lt50bp.csv")
suspicious_indels_20bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_SuspiciousIndels.lt50bp.csv")
suspicious_indels_30bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_SuspiciousIndels.lt50bp.csv")
suspicious_indels_40bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_SuspiciousIndels.lt50bp.csv")

#df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_lt50bp.csv")
#df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_lt50bp.csv")
df_3 = fread("rAF_project_2/UK.BB.exomes.430k.sites_indelsonly_rAF_lt50bp.csv")

## split varID 
# suspicious_indels_10bp_df_1 = suspicious_indels_10bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# suspicious_indels_20bp_df_1 = suspicious_indels_20bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# suspicious_indels_30bp_df_1 = suspicious_indels_30bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# suspicious_indels_40bp_df_1 = suspicious_indels_40bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# 
# suspicious_indels_10bp_df_2 = suspicious_indels_10bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# suspicious_indels_20bp_df_2 = suspicious_indels_20bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# suspicious_indels_30bp_df_2 = suspicious_indels_30bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# suspicious_indels_40bp_df_2 = suspicious_indels_40bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# 

suspicious_indels_10bp_df_3 = suspicious_indels_10bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
suspicious_indels_20bp_df_3 = suspicious_indels_20bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
suspicious_indels_30bp_df_3 = suspicious_indels_30bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
suspicious_indels_40bp_df_3 = suspicious_indels_40bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# df_1 = df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
# df_2 = df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
df_3 = df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))

# set bp ranges
bp_range = c("10", "20", "30", "40")

#df_name_1="gnomAD"
#df_name_2="IGM"
df_name_3="UK.BB"

# chromosomes to analyze 
chroms = c("X", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "Y")

## make chroms an ordered factor 
chroms <- factor (chroms, levels = chroms)

rare = 10^-4 

###########################################################################################
###########################################################################################
###########################################################################################

for (i in bp_range) {
  ## get sus for that bp range data frame 
  # df_1_sus = paste0("suspicious_indels_", i, "bp_df_1")
  # df_2_sus = paste0("suspicious_indels_", i, "bp_df_2")
  df_3_sus = paste0("suspicious_indels_", i, "bp_df_3")
  
  
  ## filter for suspicious indels within each chromosome 
  for (j in chroms) {
    df_1_sus_chrom = paste0("df_1_sus_chrom_", j, "_bp_", i)
    df_2_sus_chrom = paste0("df_2_sus_chrom_", j, "_bp_", i)
    df_3_sus_chrom = paste0("df_3_sus_chrom_", j, "_bp_", i)
    
    # assign(df_1_sus_chrom, nrow(get(df_1_sus) %>% filter(CHR == j)))
    # assign(df_2_sus_chrom, nrow(get(df_2_sus) %>% filter(CHR == j)))
    assign(df_3_sus_chrom, nrow(get(df_3_sus) %>% filter(CHR == j)))
  
  ## filter for rare sAF indels within each chromosome 
    # df_1_rare_chrom = paste0("df_1_rare_chrom_", j, "_bp_", i)
    # df_2_rare_chrom = paste0("df_2_rare_chrom_", j, "_bp_", i)
    df_3_rare_chrom = paste0("df_3_rare_chrom_", j, "_bp_", i)
    
    # assign(df_1_rare_chrom, nrow(df_1 %>% filter(sAF <= rare) %>% filter(CHR == j)))
    # assign(df_2_rare_chrom, nrow(df_2 %>% filter(sAF <= rare) %>% filter(CHR == j)))
    assign(df_3_rare_chrom, nrow(df_3 %>% filter(sAF <= rare) %>% filter(CHR == j)))
    
  ## find percentage of sus/rare * 100
    # df_1_prct = paste0("df_1_prct_chrom", j, "_bp_", i)
    # df_2_prct = paste0("df_2_prct_chrom", j, "_bp_", i)
    df_3_prct = paste0("df_3_prct_chrom", j, "_bp_", i)
    
    # assign(df_1_prct, get(df_1_sus_chrom)/get(df_1_rare_chrom) * 100 )
    # assign(df_2_prct, get(df_2_sus_chrom)/get(df_2_rare_chrom) * 100 )
    assign(df_3_prct, get(df_3_sus_chrom)/get(df_3_rare_chrom) * 100 )
    
    
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

# sus_chr_distr_df1 = data.frame(x = chroms, y1 = df_1_list_bp10, y2 = df_1_list_bp20, y3 = df_1_list_bp30, y4 = df_1_list_bp40)
# 
# sus_chr_distr_df2 = data.frame(x = chroms, y1 = df_2_list_bp10, y2 = df_2_list_bp20, y3 = df_2_list_bp30, y4 = df_2_list_bp40)

sus_chr_distr_df3 = data.frame(x = chroms, y1 = df_3_list_bp10, y2 = df_3_list_bp20, y3 = df_3_list_bp30, y4 = df_3_list_bp40)

# sus_chr_distr_plot1 = ggplot(sus_chr_distr_df1, aes(x)) +
#   geom_line(aes(y = y1, color = "10bp window"), group = 1) +
#   geom_line(aes(y = y2, color = "20bp window"), group = 1) +
#   geom_line(aes(y = y3, color = "30bp window"), group = 1) +
#   geom_line(aes(y = y4, color = "40bp window"), group = 1) + 
#   ylim (0, 30) +
#   labs (title = str_wrap(paste0("Percentage of Rare Indels that are Suspicious Across Chromosomes in the ", df_name_1, " dataset"), width = 60),
#         x = "Chromosomes",
#         y = "Percentage") 
# 
# sus_chr_distr_plot2 = ggplot(sus_chr_distr_df2, aes(x)) +
#   geom_line(aes(y = y1, color = "10bp window"), group = 1) +
#   geom_line(aes(y = y2, color = "20bp window"), group = 1) +
#   geom_line(aes(y = y3, color = "30bp window"), group = 1) +
#   geom_line(aes(y = y4, color = "40bp window"), group = 1) +
#   ylim (0,30) +
#   labs (title = str_wrap(paste0("Percentage of Rare Indels that are Suspicious Across Chromosomes in the ", df_name_2, " dataset"), width = 60),
#         x = "Chromosomes",
#         y = "Percentage")

sus_chr_distr_plot3 = ggplot(sus_chr_distr_df3, aes(x)) +
  geom_line(aes(y = y1, color = "10bp window"), group = 1) +
  geom_line(aes(y = y2, color = "20bp window"), group = 1) +
  geom_line(aes(y = y3, color = "30bp window"), group = 1) +
  geom_line(aes(y = y4, color = "40bp window"), group = 1) +
  ylim (0,30) +
  labs (title = str_wrap(paste0("Percentage of Rare Indels that are Suspicious Across Chromosomes in the ", df_name_3, " dataset"), width = 60),
        x = "Chromosomes",
        y = "Percentage")

# sus_chr_distr_plot1
# sus_chr_distr_plot2
sus_chr_distr_plot3

#summary_graph = ggarrange(sus_chr_distr_plot1, sus_chr_distr_plot2, sus_chr_distr_plot3, labels = c("A", "B", "C"), ncol = 1, nrow = 2)
summary_graph = ggarrange(sus_chr_distr_plot3, labels = c("A"), ncol = 1, nrow = 2)


summary_graph


ggsave("FigS1.jpg", width = 80, height = 60, units = c("cm"), dpi = 300)


