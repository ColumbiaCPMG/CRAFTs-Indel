# Load in packages
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
install.packages("ggbreak")
library(ggbreak)
library(tidyr)
install.packages("ggpubr")
library(ggpubr)
install.packages("gt")
library(gt)
library(gridExtra)
library(patchwork)
library

rm(list=ls())
setwd("/Users/User_1/Desktop/rAF_project_2")

db_name_1 = "gnomAD"
db_name_2 = "IGM"
db_name_3 = "UK.BB"

df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_lt50bp.csv")
df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_lt50bp.csv")
df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_lt50bp.csv")

rAF_hi_indels_10bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_20bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_30bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_40bp_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

rAF_hi_indels_10bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_20bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_30bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_40bp_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

rAF_hi_indels_10bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_20bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_30bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_indels_40bp_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

bp_range = c("10", "20", "30", "40")

df_1 = separate(df_1, col = VarID, into = c("CHR", "POS", "REF", "ALT"), sep = "-")
df_2 = separate(df_2, col = VarID, into = c("CHR", "POS", "REF", "ALT"), sep = "-")
df_3 = separate(df_3, col = VarID, into = c("CHR", "POS", "REF", "ALT"), sep = "-")

for (i in bp_range) {
  print(i)
  
  ## Find the regions that the rAF_hi indels are in
  temp_regions_1 = distinct(get(paste0("rAF_hi_indels_", i, "bp_df_1")) %>% select(paste0("GID_bp", i)))
  temp_regions_2 = distinct(get(paste0("rAF_hi_indels_", i, "bp_df_2")) %>% select(paste0("GID_bp", i)))
  temp_regions_3 = distinct(get(paste0("rAF_hi_indels_", i, "bp_df_3")) %>% select(paste0("GID_bp", i)))
  
  temp_var = paste0("GID_bp", i)
  
  temp_df_1 = paste0("filtered_df_1_", i)
  temp_df_2 = paste0("filtered_df_2_", i)
  temp_df_3 = paste0("filtered_df_3_", i)
  
  ## Keep the regions that have at least one rAF_hi indel and have to have at least 2 indels (the second indel doesn't have to be rAF_hi)
  assign(temp_df_1,  df_1 %>% filter(df_1[[temp_var]] %in% temp_regions_1[[temp_var]]) %>% group_by(get(temp_var)) %>% filter(n()>1))
  assign(temp_df_2,  df_2 %>% filter(df_2[[temp_var]] %in% temp_regions_2[[temp_var]]) %>% group_by(get(temp_var)) %>% filter(n()>1))
  assign(temp_df_3,  df_3 %>% filter(df_3[[temp_var]] %in% temp_regions_3[[temp_var]]) %>% group_by(get(temp_var)) %>% filter(n()>1))
  
  #find start and end indels 
  assign(temp_df_1, get(temp_df_1) %>% group_by(get(temp_var)) %>% mutate(start_indel = as.numeric(min(POS)), end_indel = as.numeric(max(POS))) %>% ungroup())
  assign(temp_df_2, get(temp_df_2) %>% group_by(get(temp_var)) %>% mutate(start_indel = as.numeric(min(POS)), end_indel = as.numeric(max(POS))) %>% ungroup())
  assign(temp_df_3, get(temp_df_3) %>% group_by(get(temp_var)) %>% mutate(start_indel = as.numeric(min(POS)), end_indel = as.numeric(max(POS))) %>% ungroup())
  
  ## calculate the length of region 
  assign(temp_df_1, get(temp_df_1) %>% mutate (region_length = get(temp_df_1)$end_indel - get(temp_df_1)$start_indel))
  assign(temp_df_2, get(temp_df_2) %>% mutate (region_length = get(temp_df_2)$end_indel - get(temp_df_2)$start_indel))
  assign(temp_df_3, get(temp_df_3) %>% mutate (region_length = get(temp_df_3)$end_indel - get(temp_df_3)$start_indel))
  
  ## Get a chart with the number of regions with that region length
  temp_chart_1 = paste0("region_len_1_", i)
  temp_chart_2 = paste0("region_len_2_", i)
  temp_chart_3 = paste0("region_len_3_", i)
  
  assign(temp_chart_1, distinct(get(temp_df_1) %>% select (paste0("GID_bp", i), "region_length")))
  assign(temp_chart_2, distinct(get(temp_df_2) %>% select (paste0("GID_bp", i), "region_length")))
  assign(temp_chart_3, distinct(get(temp_df_3) %>% select (paste0("GID_bp", i), "region_length")))
  
  ## Get a graph 
  graph_df1 = paste0("graph_df1_bp", i)
  graph_df2 = paste0("graph_df2_bp", i)
  graph_df3 = paste0("graph_df3_bp", i)
  
  
  ## Make graphs 
  if (i == "10"){
    assign(graph_df1, ggplot(get(temp_chart_1), aes(x = region_length)) +
            geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
            geom_density(fill="grey", alpha = .5) +
            scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) +
            scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 900, 1000)) +
            ylim(0, 0.3) +
            labs(title = paste0(i, "bp region"), x = "Region Lengths (bps)", y = "Density") +
            theme(axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text=element_text(size=20), axis.text.x = element_text(angle = 90), axis.title=element_text(size=20), plot.title=element_text(size=32, hjust = 0.5, margin = margin( 1, 1, 1, 1, "cm")), plot.margin = margin( 2, 2, 2, 2, "cm"), axis.title.y = element_text(margin = margin( 1, 1, 1, 1, "cm"), angle = 90)))
  
    assign(graph_df2, ggplot(get(temp_chart_2), aes(x = region_length)) +
             geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
             geom_density(fill="grey", alpha = .5) +
             scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) +
             scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 900, 1000)) +
             ylim(0, 0.3) +
             labs(title = paste0(" "), x = "Region Lengths (bps)", y = "Density") +
             theme(axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text=element_text(size=20),  axis.text.x = element_text(angle = 90), axis.title=element_text(size=20), plot.title=element_text(size=32, hjust = 0.5, margin = margin( 1, 1, 1, 1, "cm")), plot.margin = margin( 2, 2, 2, 2, "cm"), axis.title.y = element_text(margin = margin( 1, 1, 1, 1, "cm"), angle = 90)))
  
    assign(graph_df3, ggplot(get(temp_chart_3), aes(x = region_length)) +
            geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
            geom_density(fill="grey", alpha = .5) + 
            scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) + 
            scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 900, 1000)) +
            ylim(0, 0.3) + 
            labs(title = paste0(" "), x = "Region Lengths (bps)", y = "Density") + 
            theme(axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text=element_text(size=20),  axis.text.x = element_text(angle = 90), axis.title=element_text(size=20), plot.title=element_text(size=32, hjust = 0.5, margin = margin( 1, 1, 1, 1, "cm")), plot.margin = margin( 2, 2, 2, 2, "cm"), axis.title.y = element_text(margin = margin( 1, 1, 1, 1, "cm"), angle = 90)))
  }
  #makes graphs without y axis
  if ( i != "10"){
    assign(graph_df1, ggplot(get(temp_chart_1), aes(x = region_length)) +
             geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
             geom_density(fill="grey", alpha = .5) +
             scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) +
             scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 900, 1000)) +
             ylim(0, 0.3) +
             labs(title = paste0(i, "bp region"), x = "Region Lengths (bps)") +
             theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text=element_text(size=20), axis.text.x = element_text(angle = 90), axis.title=element_text(size=20), plot.title=element_text(size=32, hjust = 0.5, margin = margin( 1, 1, 1, 1, "cm")), plot.margin = margin( 2, 2, 2, 2, "cm")))
    assign(graph_df2, ggplot(get(temp_chart_2), aes(x = region_length)) +
             geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
             geom_density(fill="grey", alpha = .5) +
             scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) +
             scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 900, 1000)) +
             ylim(0, 0.3) +
             labs(title = paste0(" "), x = "Region Lengths (bps)") +
             theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text=element_text(size=20),  axis.text.x = element_text(angle = 90), axis.title=element_text(size=20), plot.title=element_text(size=32, hjust = 0.5, margin = margin( 1, 1, 1, 1, "cm")), plot.margin = margin( 2, 2, 2, 2, "cm")))
    
    assign(graph_df3, ggplot(get(temp_chart_3), aes(x = region_length)) +
             geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
             geom_density(fill="grey", alpha = .5) + 
             scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) + 
             scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 900, 1000)) +
             ylim(0, 0.3) + 
             labs(title = paste0(" "), x = "Region Lengths (bps)") + 
             theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),axis.text.x.top = element_blank(),axis.ticks.x.top = element_blank(), axis.line.x.top = element_blank(), axis.text=element_text(size=20),  axis.text.x = element_text(angle = 90), axis.title=element_text(size=20), plot.title=element_text(size=32, hjust = 0.5, margin = margin( 1, 1, 1, 1, "cm")), plot.margin = margin( 2, 2, 2, 2, "cm")))
    
    
    
    
  }
}

# manually scales some of the graphs beacuse the graphs with the y axis would appear visually smaller then the graphs without the y axis
scaled_graph_df1_bp10 <- print(graph_df1_bp10) + theme(plot.margin = margin(0, .2,  0, -2.7, "cm"))
scaled_graph_df1_bp20 <- print(graph_df1_bp20) + theme(plot.margin = margin( 0, .2, 0, 0, "cm"))
scaled_graph_df1_bp30 <- print(graph_df1_bp30) + theme(plot.margin = margin( 0, .2, 0, 0, "cm"))
scaled_graph_df1_bp40 <- print(graph_df1_bp40) + theme(plot.margin = margin( 0, .2, 0, 0, "cm"))

scaled_graph_df2_bp10 <- print(graph_df2_bp10) + theme(plot.margin = margin(0, .2,  0, -2.7, "cm"))
scaled_graph_df2_bp20 <- print(graph_df2_bp20) + theme(plot.margin = margin( 0, .2, 0, 0, "cm"))
scaled_graph_df2_bp30 <- print(graph_df2_bp30) + theme(plot.margin = margin( 0, .2, 0, 0, "cm"))
scaled_graph_df2_bp40 <- print(graph_df2_bp40) + theme(plot.margin = margin( 0, .2, 0, 0, "cm"))

scaled_graph_df3_bp10 <- print(graph_df3_bp10) + theme(plot.margin = margin(0, .2,  0, -2.7, "cm"))
scaled_graph_df3_bp20 <- print(graph_df3_bp20) + theme(plot.margin = margin( 0, .2, 0, 0, "cm"))
scaled_graph_df3_bp30 <- print(graph_df3_bp30) + theme(plot.margin = margin( 0, .2, 0, 0, "cm"))
scaled_graph_df3_bp40 <- print(graph_df3_bp40) + theme(plot.margin = margin( 0, .2, 0, 0, "cm"))




figure_ = ggarrange(scaled_graph_df1_bp10,
          scaled_graph_df1_bp20,
          scaled_graph_df1_bp30,
          scaled_graph_df1_bp40,
          scaled_graph_df2_bp10,
          scaled_graph_df2_bp20,
          scaled_graph_df2_bp30,
          scaled_graph_df2_bp40,
          scaled_graph_df3_bp10,
          scaled_graph_df3_bp20,
          scaled_graph_df3_bp30,
          scaled_graph_df3_bp40,
          labels = c(" ", " ", " ", " ", " IGM", " ", " ", " ", "UK.BB", " ", " ", " "),
          ncol = 4, nrow = 3, font.label = list(size = 32, color = "black")) + theme(plot.margin = unit(c(2,2,2,2), "cm"))



# added the gnomAD label this way so it would not overlap with the grpah title
labled_figure <- annotate_figure(figure_,
                               top = text_grob("gnomAD", 
                                               color = "black", 
                                               face = "bold", 
                                               size = 32, 
                                               hjust = 7.3,
                                               vjust = 4))

# Adjust the margins of the plot for better appearance
labled_figure + theme(plot.margin = unit(c(2, 2, 2, 2), "cm"))
ggsave("Fig2.jpg", width = 80, height = 60, units = c("cm"), dpi = 300)



## To find the mean and median region lengths.
for (i in bp_range) {
  # get the dataframe with all the indels in a region with a rAF_hi indel
  temp_var1 = paste0("filtered_df_1_", i)
  temp_var2 = paste0("filtered_df_2_", i)
  temp_var3 = paste0("filtered_df_3_", i)

  # find mean and median
  mean_df1 = paste0("mean_region_len_df1_bp", i)
  mean_df2 = paste0("mean_region_len_df2_bp", i)
  mean_df3 = paste0("mean_region_len_df3_bp", i)

  gID = paste0("GID_bp", i)


  ## remember: got to collapse by gID
  assign(mean_df1, mean(unique(get(temp_var1) %>% select(gID, region_length))$region_length))
  assign(mean_df2, mean(unique(get(temp_var2) %>% select(gID, region_length))$region_length))
  assign(mean_df3, mean(unique(get(temp_var3) %>% select(gID, region_length))$region_length))

  median_df1 = paste0("median_region_len_df1_bp", i)
  median_df2 = paste0("median_region_len_df2_bp", i)
  median_df3 = paste0("median_region_len_df3_bp", i)

  assign(median_df1, median(unique(get(temp_var1) %>% select(gID, region_length))$region_length))
  assign(median_df2, median(unique(get(temp_var2) %>% select(gID, region_length))$region_length))
  assign(median_df3, median(unique(get(temp_var3) %>% select(gID, region_length))$region_length))

  print(paste0("Mean region length for ", db_name_1, " and bp region ", i, ": ", get(mean_df1)))
  print(paste0("Mean region length for ", db_name_2, " and bp region ", i, ": ", get(mean_df2)))
  print(paste0("Mean region length for ", db_name_3, " and bp region ", i, ": ", get(mean_df3)))

  print(paste0("Median region length for ", db_name_1, " and bp region ", i, ": ", get(median_df1)))
  print(paste0("Median region length for ", db_name_2, " and bp region ", i, ": ", get(median_df2)))
  print(paste0("Median region length for ", db_name_3, " and bp region ", i, ": ", get(median_df3)))
}



for (i in bp_range) {

  ## get the region_len dataframe
  ben_len_1 = paste0("region_len_1_", i)
  ben_len_2 = paste0("region_len_2_", i)
  ben_len_3 = paste0("region_len_3_", i)

  ## assign the count region variable
  count_region_1 = paste0("count_region", i, "_df1")
  count_region_2 = paste0("count_region", i, "_df2")
  count_region_3 = paste0("count_region", i, "_df3")

  assign(count_region_1, get(ben_len_1) %>% group_by(region_length) %>% count(region_length))
  assign(count_region_2, get(ben_len_2) %>% group_by(region_length) %>% count(region_length))
  assign(count_region_3, get(ben_len_3) %>% group_by(region_length) %>% count(region_length))

  ## find the percentage
  percent_df1 = paste0("percent_lt_", i, "bp_region_", i, "df1")
  percent_df2 = paste0("percent_lt_", i, "bp_region_", i, "df2")
  percent_df3 = paste0("percent_lt_", i, "bp_region_", i, "df3")

  assign(percent_df1, sum((get(count_region_1) %>% filter(region_length < as.numeric(i)))$n) / sum(get(count_region_1)$n) * 100 )
  assign(percent_df2, sum((get(count_region_2) %>% filter(region_length < as.numeric(i)))$n) / sum(get(count_region_2)$n) * 100 )
  assign(percent_df3, sum((get(count_region_3) %>% filter(region_length < as.numeric(i)))$n) / sum(get(count_region_3)$n) * 100 )

  print(paste0("Percent less than ", i, " bps, region ", i, " is ", get(percent_df1), "% for dataset ", db_name_1))
  print(paste0("Percent less than ", i, " bps, region ", i, " is ", get(percent_df2), "% for dataset ", db_name_2))
  print(paste0("Percent less than ", i, " bps, region ", i, " is ", get(percent_df3), "% for dataset ", db_name_3))

}



for (i in bp_range) {
  # get the dataframe with all the indels in a region with a rAF_hi indel
  temp_var1 = paste0("filtered_df_1_", i)
  temp_var2 = paste0("filtered_df_2_", i)
  temp_var3 = paste0("filtered_df_3_", i)

  # assign a variable for the regions
  regions_df1 = paste0 ("regions_", i, "df1")
  regions_df2 = paste0 ("regions_", i, "df2")
  regions_df3 = paste0 ("regions_", i, "df3")

  assign (regions_df1, nrow(unique(get(temp_var1)[paste0("GID_bp", i)])))
  assign (regions_df2, nrow(unique(get(temp_var2)[paste0("GID_bp", i)])))
  assign (regions_df3, nrow(unique(get(temp_var3)[paste0("GID_bp", i)])))

  print(paste0("Number of regions with rAF_hi indels in ", db_name_1, " is: ", get(regions_df1)))
  print(paste0("Number of regions with rAF_hi indels in ", db_name_2, " is: ", get(regions_df2)))
  print(paste0("Number of regions with rAF_hi indels in ", db_name_3, " is: ", get(regions_df3)))
}


df_1_regions_rAF_hi_indels = c(regions_10df1, regions_20df1, regions_30df1, regions_40df1)
df_1_mean_region_len = c(mean_region_len_df1_bp10, mean_region_len_df1_bp20, mean_region_len_df1_bp30, mean_region_len_df1_bp40)
df_1_median_region_len = c(median_region_len_df1_bp10, median_region_len_df1_bp20, median_region_len_df1_bp30, median_region_len_df1_bp40)
df_1_regions_prct_under_region_size = c(percent_lt_10bp_region_10df1, percent_lt_20bp_region_20df1, percent_lt_30bp_region_30df1, percent_lt_40bp_region_40df1)

df_1_summary = data.frame(bp_range, df_1_regions_rAF_hi_indels, df_1_mean_region_len, df_1_median_region_len, df_1_regions_prct_under_region_size)

#rounds numbers in the data frame
df_1_summary <- df_1_summary %>%
  mutate_at(vars(df_1_regions_rAF_hi_indels, df_1_mean_region_len,df_1_regions_prct_under_region_size), ~round(., 2))


# ads percent sign to the percent column
df_1_summary <- df_1_summary %>%
  mutate_at(vars(df_1_regions_prct_under_region_size), ~paste0(., "%"))

#adds commas to the numbers in the data frame
df_1_summary <- df_1_summary %>%
  mutate(df_1_regions_rAF_hi_indels = format(df_1_regions_rAF_hi_indels, big.mark = ","))

df_1_summary <- df_1_summary %>%
  mutate(bp_range = paste(bp_range = paste(bp_range,"bps")))

#renames the columns in teh data frame
colnames(df_1_summary) <- c("Max region length", "Nb. of rAF-hi genomic regions", "Mean length of rAF-hi genomic region (bp)", "Median length of rAF-hi suspicious genomic region (bp)","Proportion of rAF-hi genomic regions smaller than the max region length")
#df_1_summary <- knitr::kable(df_1_summary, caption = "<b>Table S3:</b> Length of suspicious regions in the gnomAD dataset")


# sets a limit on the column rows, so that the titles dont make the columns wider
#this appears in 2 other places
# max_width_per_line <- 24
#
#
# col_names <- names(df_1_summary)
# col_names <- lapply(col_names, function(x) {
#   split_lines <- strwrap(x, width = max_width_per_line)
#   paste(split_lines, collapse = "\n")
# })

#names(df_1_summary) <- col_names
#options(width = 20)

# df_1_summary <- df_1_summary %>%
#   tab_header(title = md("Data listing from **gtcars**"))

df_1_summary



#####df2 ##########

df_2_regions_rAF_hi_indels = c(regions_10df2, regions_20df2, regions_30df2, regions_40df2)
df_2_mean_region_len = c(mean_region_len_df2_bp10, mean_region_len_df2_bp20, mean_region_len_df2_bp30, mean_region_len_df2_bp40)
df_2_median_region_len = c(median_region_len_df2_bp10, median_region_len_df2_bp20, median_region_len_df2_bp30, median_region_len_df2_bp40)
df_2_regions_prct_under_region_size = c(percent_lt_10bp_region_10df2, percent_lt_20bp_region_20df2, percent_lt_30bp_region_30df2, percent_lt_40bp_region_40df2)

df_2_summary = data.frame(bp_range, df_2_regions_rAF_hi_indels, df_2_mean_region_len, df_2_median_region_len, df_2_regions_prct_under_region_size)
#rounds numbers
df_2_summary <- df_2_summary %>%
  mutate_at(vars(df_2_regions_rAF_hi_indels, df_2_mean_region_len,df_2_regions_prct_under_region_size), ~round(., 2))
# ads percent sign to the percent column
df_2_summary <- df_2_summary %>%
  mutate_at(vars(df_2_regions_prct_under_region_size), ~paste0(., "%"))
df_2_summary <- df_2_summary %>%
  mutate(df_2_regions_rAF_hi_indels = format(df_2_regions_rAF_hi_indels, big.mark = ","))

df_2_summary <- df_2_summary %>%
  mutate(bp_range = paste(bp_range = paste(bp_range,"bps")))

colnames(df_2_summary) <- c("Max region length", "Nb. of rAF-hi genomic regions", "Mean length of rAF-hi genomic region (bp)", "Median length of rAF-hi suspicious genomic region (bp)","Proportion of rAF-hi genomic regions smaller than the max region length")

# sets a limit on the column rows, so that the titles dont make the columns wider

# max_width_per_line <- 20
#
#
# col_names <- names(df_2_summary)
# col_names <- lapply(col_names, function(x) {
#   split_lines <- strwrap(x, width = max_width_per_line)
#   paste(split_lines, collapse = "\n")
# })
#
# names(df_2_summary) <- col_names
df_2_summary


####df3#######


df_3_regions_rAF_hi_indels = c(regions_10df3, regions_20df3, regions_30df3, regions_40df3)
df_3_mean_region_len = c(mean_region_len_df3_bp10, mean_region_len_df3_bp20, mean_region_len_df3_bp30, mean_region_len_df3_bp40)
df_3_median_region_len = c(median_region_len_df3_bp10, median_region_len_df3_bp20, median_region_len_df3_bp30, median_region_len_df3_bp40)
df_3_regions_prct_under_region_size = c(percent_lt_10bp_region_10df3, percent_lt_20bp_region_20df3, percent_lt_30bp_region_30df3, percent_lt_40bp_region_40df3)

df_3_summary = data.frame(bp_range, df_3_regions_rAF_hi_indels, df_3_mean_region_len, df_3_median_region_len, df_3_regions_prct_under_region_size)
df_3_summary <- df_3_summary %>%
  mutate_at(vars(df_3_regions_rAF_hi_indels, df_3_mean_region_len,df_3_regions_prct_under_region_size), ~round(., 2))
# ads percent sign to the percent column
df_3_summary <- df_3_summary %>%
  mutate_at(vars(df_3_regions_prct_under_region_size), ~paste0(., "%"))

df_3_summary <- df_3_summary %>%
  mutate(df_3_regions_rAF_hi_indels = format(df_3_regions_rAF_hi_indels, big.mark = ","))
# Create the table with the title

df_3_summary <- df_3_summary %>%
  mutate(bp_range = paste(bp_range = paste(bp_range,"bps")))
colnames(df_3_summary) <- c("Range", "Nb. of rAF-hi genomic regions", "Mean length of rAF-hi genomic region (bp)", "Median length of rAF-hi genomic region (bp)","Proportion of rAF-hi genomic regions smaller than the range used")


# sets a limit on the column rows, so that the titles dont make the columns wider

# max_width_per_line <- 20
#
#
# col_names <- names(df_3_summary)
# col_names <- lapply(col_names, function(x) {
#   split_lines <- strwrap(x, width = max_width_per_line)
#   paste(split_lines, collapse = "\n")
# })
#
# names(df_3_summary) <- col_names
df_3_summary



fwrite(df_1_summary, "Table S3: Length of suspicious regions in the gnomAD dataset.csv")
fwrite(df_2_summary, "Table S4: Length of suspicious regions in the IGM dataset.csv")
fwrite(df_3_summary, "Table S5: Length of suspicious regions in the UK.BB dataset.csv")



###############################################
###############################################
###############################################

# P- value calculations for tables S4



#xx$region_length
dfs = c("df_1_", "df_2_", "df_3_")

for (i in dfs){
    df_10_bp = get(paste0("filtered_",i,"10"))
    df_40_bp = get(paste0("filtered_",i,"40"))
    
    df_10_bp_uniques = paste0("df_10_bp_uniques_", i)
    df_40_bp_uniques = paste0("df_40_bp_uniques_", i)
    
    assign( df_10_bp_uniques, (unique(df_10_bp %>% select(GID_bp10, region_length))$region_length))
    assign(df_40_bp_uniques,  (unique(df_40_bp %>% select(GID_bp40, region_length))$region_length))
    
    
   # bp10_lengths = df_10_bp_uniques$region_length
  #  bp40_lengths = df_40_bp_uniques$region_length
    
    print(i)
    options(digits = 2)
    print(paste0("Mean 10bp is", median(get(df_10_bp_uniques))))
    print(paste0("Mean 40bp is", median(get(df_40_bp_uniques))))
    result = t.test(get(df_10_bp_uniques), get(df_40_bp_uniques), var.equal = TRUE )
    print(result)
   # print(result, digits = 4)
    print("_______________")
  
}  

#df_export = cbind(as.data.frame(df_10_bp_uniques_df_1_), as.data.frame(df_40_bp_uniques_df_1_))

merged_df <- merge(as.data.frame(df_10_bp_uniques_df_1_), as.data.frame(df_40_bp_uniques_df_1_), all = TRUE)

fwrite(region_len_1_10, "gnomad_10bp_region_lengths.csv")
fwrite(region_len_1_40, "gnomad_40bp_region_lengths.csv")
