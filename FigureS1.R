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
install.packages("ggh4x")
library(ggh4x)


###########################################################################################
###########################################################################################
###########################################################################################

setwd("")

## Read in files 
rAF_hi_indels_10bp_df_1 = fread("gnomAD/gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_20bp_df_1 = fread("gnomAD/gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_30bp_df_1 = fread("gnomAD/gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_40bp_df_1 = fread("gnomAD/gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

rAF_hi_indels_10bp_df_2 = fread("IGM/2023-03-23_IGM_n39367_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_20bp_df_2 = fread("IGM/2023-03-23_IGM_n39367_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_30bp_df_2 = fread("IGM/2023-03-23_IGM_n39367_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_40bp_df_2 = fread("IGM/2023-03-23_IGM_n39367_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

rAF_hi_indels_10bp_df_3 = fread("UKBB/UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_20bp_df_3 = fread("UKBB/UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_30bp_df_3 = fread("UKBB/UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_40bp_df_3 = fread("UKBB/UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

df_1 = fread("gnomAD/gnomad.exomes.r2.1.1.sites_indelsonly_rAF_lt50bp.csv")
df_2 = fread("IGM/2023-03-23_IGM_n39367_indelsonly_rAF_lt50bp.csv")
df_3 = fread("UKBB/UK.BB.exomes.430k.sites_indelsonly_rAF_lt50bp.csv")

## split varID 
rAF_hi_indels_10bp_df_1 = rAF_hi_indels_10bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_20bp_df_1 = rAF_hi_indels_20bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_30bp_df_1 = rAF_hi_indels_30bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_40bp_df_1 = rAF_hi_indels_40bp_df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
#
rAF_hi_indels_10bp_df_2 = rAF_hi_indels_10bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_20bp_df_2 = rAF_hi_indels_20bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_30bp_df_2 = rAF_hi_indels_30bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_40bp_df_2 = rAF_hi_indels_40bp_df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
#

rAF_hi_indels_10bp_df_3 = rAF_hi_indels_10bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_20bp_df_3 = rAF_hi_indels_20bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_30bp_df_3 = rAF_hi_indels_30bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
rAF_hi_indels_40bp_df_3 = rAF_hi_indels_40bp_df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
df_1 = df_1 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
df_2 = df_2 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))
df_3 = df_3 %>% separate(VarID, c("CHR", "POS", "REF", "ALT"))

# set bp ranges
bp_range = c("10", "20", "30", "40")

df_name_1="gnomAD"
df_name_2="IGM"
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
  df_1_rAF_hi = paste0("rAF_hi_indels_", i, "bp_df_1")
  df_2_rAF_hi = paste0("rAF_hi_indels_", i, "bp_df_2")
  df_3_rAF_hi = paste0("rAF_hi_indels_", i, "bp_df_3")
  
  
  ## filter for rAF_hi indels within each chromosome 
  for (j in chroms) {
    df_1_rAF_hi_chrom = paste0("df_1_rAF_hi_chrom_", j, "_bp_", i)
    df_2_rAF_hi_chrom = paste0("df_2_rAF_hi_chrom_", j, "_bp_", i)
    df_3_rAF_hi_chrom = paste0("df_3_rAF_hi_chrom_", j, "_bp_", i)
    
    assign(df_1_rAF_hi_chrom, nrow(get(df_1_rAF_hi) %>% filter(CHR == j)))
    assign(df_2_rAF_hi_chrom, nrow(get(df_2_rAF_hi) %>% filter(CHR == j)))
    assign(df_3_rAF_hi_chrom, nrow(get(df_3_rAF_hi) %>% filter(CHR == j)))
    
    ## filter for rAF_lo sAF indels within each chromosome 
    df_1_rAF_lo_chrom = paste0("df_1_rAF_lo_chrom_", j, "_bp_", i)
    df_2_rAF_lo_chrom = paste0("df_2_rAF_lo_chrom_", j, "_bp_", i)
    df_3_rAF_lo_chrom = paste0("df_3_rAF_lo_chrom_", j, "_bp_", i)
    
    assign(df_1_rAF_lo_chrom, nrow(df_1 %>% filter(sAF <= rAF_lo) %>% filter(CHR == j)))
    assign(df_2_rAF_lo_chrom, nrow(df_2 %>% filter(sAF <= rAF_lo) %>% filter(CHR == j)))
    assign(df_3_rAF_lo_chrom, nrow(df_3 %>% filter(sAF <= rAF_lo) %>% filter(CHR == j)))
    
    ## find percentage of rAF_hi/rAF_lo * 100
    df_1_prct = paste0("df_1_prct_chrom", j, "_bp_", i)
    df_2_prct = paste0("df_2_prct_chrom", j, "_bp_", i)
    df_3_prct = paste0("df_3_prct_chrom", j, "_bp_", i)
    
    assign(df_1_prct, get(df_1_rAF_hi_chrom)/get(df_1_rAF_lo_chrom) * 100 )
    assign(df_2_prct, get(df_2_rAF_hi_chrom)/get(df_2_rAF_lo_chrom) * 100 )
    assign(df_3_prct, get(df_3_rAF_hi_chrom)/get(df_3_rAF_lo_chrom) * 100 )
    
    
  }
  
  ## make a list for graph data frame 
  df_1_list = paste0("df_1_list_bp", i)
  df_2_list = paste0("df_2_list_bp", i)
  df_3_list = paste0("df_3_list_bp", i)
  
  assign(df_1_list, c())
  assign(df_2_list, c())
  assign(df_3_list, c())
  
  for (k in chroms) {
    ## get variable 
    get_temp_df1 = paste0("df_1_prct_chrom", k, "_bp_", i)
    get_temp_df2 = paste0("df_2_prct_chrom", k, "_bp_", i)
    get_temp_df3 = paste0("df_3_prct_chrom", k, "_bp_", i)
    
    assign(df_1_list, c(get(df_1_list), get(get_temp_df1)))
    assign(df_2_list, c(get(df_2_list), get(get_temp_df2)))
    assign(df_3_list, c(get(df_3_list), get(get_temp_df3)))
  }
}

###########################################################################################
###########################################################################################
###########################################################################################

rAF_hi_chr_distr_df1 = data.frame(x = chroms, y1 = df_1_list_bp10, y2 = df_1_list_bp20, y3 = df_1_list_bp30, y4 = df_1_list_bp40)

rAF_hi_chr_distr_df2 = data.frame(x = chroms, y1 = df_2_list_bp10, y2 = df_2_list_bp20, y3 = df_2_list_bp30, y4 = df_2_list_bp40)

rAF_hi_chr_distr_df3 = data.frame(x = chroms, y1 = df_3_list_bp10, y2 = df_3_list_bp20, y3 = df_3_list_bp30, y4 = df_3_list_bp40)


###########################################################################################
###########################################################################################
###########################################################################################

#Create the merged plot

rAF_hi_chr_distr_df1$Plot <- "gnomAD"
rAF_hi_chr_distr_df2$Plot <- "IGM"
rAF_hi_chr_distr_df3$Plot <- "UK.BB"



combined_data <- rbind( rAF_hi_chr_distr_df1, rAF_hi_chr_distr_df2, rAF_hi_chr_distr_df3)

# ensures that the legend contains the proper variable names 

combined_data_long <- combined_data %>%
  gather(key = "variable", value = "value", -x, -Plot)%>%
  mutate(variable = case_when(
    variable == "y1" ~ "10bp",
    variable == "y2" ~ "20bp",
    variable == "y3" ~ "30bp",
    variable == "y4" ~ "40bp",
    TRUE ~ variable
  ))

options(repr.plot.width = 10, repr.plot.height = 8)

df_scales <- data.frame(
  Panel = c("gnomad", "IGM", "UK.BB"), 
  ymin = c(0, 0, 0), #sets the y min for each plot
  ymax = c(30, 30, 30), #sets the y max for each plot 
  n = c(2, 2, 2) #sets the number of axis ticks displayed 
)
df_scales <- split(df_scales, df_scales$Panel) #makes the y lim scales independent for each panel 


#applys df_scales to 'scales', used when plotting combined_plot
scales <- lapply(df_scales, function(x) {
  scale_y_continuous(limits = c(x$ymin, x$ymax), n.breaks = x$n)
}) 

combined_plot<- ggplot(combined_data_long, aes(x, value, color = variable, group = variable)) +
  geom_line() +
  facet_grid(Plot ~ ., scales = "free_y", space = "free_y", switch = "y") +
  ggh4x::facetted_pos_scales(
    y = scales
  )+
labs(x = "Chromosome", y = "Percentage", title = "Percentage of sAF-lo Indels that are rAF-hi Across Chromosomes in the gnomAD, IGM, and UKBB Data Sets.") +
  labs(color = "Range") +
  #removes grid lines 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

setwd("")

ggsave("FigS1.jpg", width = 40, height = 20, units = c("cm"), dpi = 300)
