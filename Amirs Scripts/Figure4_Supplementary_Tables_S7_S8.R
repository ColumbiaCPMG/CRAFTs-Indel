library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
install.packages("ggpubr")
library(ggpubr)
install.packages("ggfittext")
library(ggfittext)
install.packages("ggrepel")
library(ggrepel)

rm(list=ls())


rAF_lo = (1 * 10^-4)

bp_range = c("10", "20", "30", "40")

conditions = c("rAF_hi", "rAF_lo", "sAF_hi")

output_path = ""
setwd("/Users/User_1/Desktop/rAF_project_2/")

df_name_1 = "gnomAD"
df_name_2 = "IGM"
df_name_3 = "UK.BB"

rAF_hi_bp_10_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_20_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_30_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_40_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.region.bed")

rAF_hi_bp_10_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_20_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_30_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_40_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.region.bed")

rAF_hi_bp_10_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_20_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_30_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_40_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.region.bed")


rAF_lo_bp_10_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_rAF_loIndels.lt50bp.region.bed")
rAF_lo_bp_20_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_rAF_loIndels.lt50bp.region.bed")
rAF_lo_bp_30_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_rAF_loIndels.lt50bp.region.bed")
rAF_lo_bp_40_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_rAF_loIndels.lt50bp.region.bed")

rAF_lo_bp_10_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_rAF_loIndels.lt50bp.region.bed")
rAF_lo_bp_20_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp20_rAF_loIndels.lt50bp.region.bed")
rAF_lo_bp_30_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp30_rAF_loIndels.lt50bp.region.bed")
rAF_lo_bp_40_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp40_rAF_loIndels.lt50bp.region.bed")

rAF_lo_bp_10_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_rAF_loIndels.lt50bp.region.bed")
rAF_lo_bp_20_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_rAF_loIndels.lt50bp.region.bed")
rAF_lo_bp_30_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_rAF_loIndels.lt50bp.region.bed")
rAF_lo_bp_40_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_rAF_loIndels.lt50bp.region.bed")

sAF_hi_bp_10_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_sAF_hiIndels.lt50bp.region.bed")
sAF_hi_bp_20_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_sAF_hiIndels.lt50bp.region.bed")
sAF_hi_bp_30_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_sAF_hiIndels.lt50bp.region.bed")
sAF_hi_bp_40_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_sAF_hiIndels.lt50bp.region.bed")

sAF_hi_bp_10_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_sAF_hiIndels.lt50bp.region.bed")
sAF_hi_bp_20_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp20_sAF_hiIndels.lt50bp.region.bed")
sAF_hi_bp_30_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp30_sAF_hiIndels.lt50bp.region.bed")
sAF_hi_bp_40_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp40_sAF_hiIndels.lt50bp.region.bed")

sAF_hi_bp_10_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_sAF_hiIndels.lt50bp.region.bed")
sAF_hi_bp_20_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_sAF_hiIndels.lt50bp.region.bed")
sAF_hi_bp_30_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_sAF_hiIndels.lt50bp.region.bed")
sAF_hi_bp_40_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_sAF_hiIndels.lt50bp.region.bed")

system("bp_range=('10' '20' '30' '40')
condition=('rAF_hi' 'rAF_lo' 'sAF_hi')

df_name_1='gnomAD'
df_name_2='IGM'
df_name_3='UK.BB'

my_date=''

output_path='/Users/User_1/Desktop/rAF_project_2/'


lcr_file_37='numInsteadOfGRCh37_AllTandemRepeatsandHomopolymers_slop5.bed'
lcr_file_38='numInsteadOfChr_GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed'

for n in ${condition[@]};
do 
   for i in ${bp_range[@]};
   do
      echo $i
      bedtools intersect -wa -a ${output_path}gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp${i}_${n}Indels.lt50bp.region.bed -b ${output_path}${lcr_file_37} | sort -k1,1n -k2,2n | uniq > ${output_path}${my_date}_LCR_intersect_${n}_${i}_df_1_${df_name_1}.bed
      bedtools intersect -wa -a ${output_path}2023-03-23_IGM_n39367_indelsonly_rAF_bp${i}_${n}Indels.lt50bp.region.bed -b ${output_path}${lcr_file_37} | sort -k1,1n -k2,2n | uniq > ${output_path}${my_date}_LCR_intersect_${n}_${i}_df_2_${df_name_2}.bed
      bedtools intersect -wa -a ${output_path}UK.BB.exomes.430k.sites_indelsonly_rAF_bp${i}_${n}Indels.lt50bp.region.bed -b ${output_path}${lcr_file_38} | sort -k1,1n -k2,2n | uniq > ${output_path}${my_date}_LCR_intersect_${n}_${i}_df_3_${df_name_3}.bed
   done;
done", intern = TRUE)

formatter10000 = function(x) {
  x/10000
}
my_date = ''
for (i in conditions) {
  for (j in bp_range) {
    # get files with the regions before overlap with LCR 
    total_file_df1 = paste0(i, "_bp_", j, "_df_1")
    total_file_df2 = paste0(i, "_bp_", j, "_df_2")
    total_file_df3 = paste0(i, "_bp_", j, "_df_3")
    
    # get files with regions after overlap with LCR 
    LCR_file_df1 = paste0(output_path, my_date, "_LCR_intersect_", i,"_", j, "_df_1_", df_name_1, ".bed")
    LCR_file_df2 = paste0(output_path, my_date, "_LCR_intersect_", i,"_", j, "_df_2_", df_name_2, ".bed")
    LCR_file_df3 = paste0(output_path, my_date, "_LCR_intersect_", i,"_", j, "_df_3_", df_name_3, ".bed")
    
    n_total_df1 = paste0(i, "_total_bp_", j, "_df_1")
    n_total_df2 = paste0(i, "_total_bp_", j, "_df_2")
    n_total_df3 = paste0(i, "_total_bp_", j, "_df_3")
    
    assign( n_total_df1 , nrow(get(total_file_df1)))
    assign( n_total_df2 , nrow(get(total_file_df2)))
    assign( n_total_df3 , nrow(get(total_file_df3)))
    
    # print(paste0(n_total_df1, ": ", get(n_total_df1)))
    # print(paste0(n_total_df2, ": ", get(n_total_df2)))
    # print(paste0(n_total_df3, ": ", get(n_total_df3)))
    
    
    n_LCR_df1 = paste0(i, "_LCR_bp_", j, "_df_1")
    n_LCR_df2 = paste0(i, "_LCR_bp_", j, "_df_2")
    n_LCR_df3 = paste0(i, "_LCR_bp_", j, "_df_3")
    
    assign( n_LCR_df1, nrow(fread(LCR_file_df1)))
    assign( n_LCR_df2, nrow(fread(LCR_file_df2)))
    assign( n_LCR_df3, nrow(fread(LCR_file_df3)))
    
    # print(paste0(n_LCR_df1, ": ", get(n_LCR_df1)))
    # print(paste0(n_LCR_df2, ": ", get(n_LCR_df2)))
    # print(paste0(n_LCR_df3, ": ", get(n_LCR_df3)))
    
    # get not in lcr 
    n_not_LCR_df1 = paste0(i, "_not_in_LCR_bp_", j,"_df_1")
    n_not_LCR_df2 = paste0(i, "_not_in_LCR_bp_", j,"_df_2")
    n_not_LCR_df3 = paste0(i, "_not_in_LCR_bp_", j,"_df_3")
    
    assign( n_not_LCR_df1, get(n_total_df1) - get(n_LCR_df1))
    assign( n_not_LCR_df2, get(n_total_df2) - get(n_LCR_df2))
    assign( n_not_LCR_df3, get(n_total_df3) - get(n_LCR_df3))
    
    #print(paste0(n_not_LCR_df1, ": ", get(n_not_LCR_df1)))
    #print(paste0(n_not_LCR_df2, ": ", get(n_not_LCR_df2)))
    
    # get percent 
    prct_LCR_df1 = paste0(i, "_prct_LCR_bp_", j, "_df_1")
    prct_LCR_df2 = paste0(i, "_prct_LCR_bp_", j, "_df_2")
    prct_LCR_df3 = paste0(i, "_prct_LCR_bp_", j, "_df_3")
    
    assign( prct_LCR_df1, paste0(round(get(n_LCR_df1)/get(n_total_df1) * 100, 0), "%"))
    assign( prct_LCR_df2, paste0(round(get(n_LCR_df2)/get(n_total_df2) * 100, 0), "%"))
    assign( prct_LCR_df3, paste0(round(get(n_LCR_df3)/get(n_total_df3) * 100, 0), "%"))
    
    print(paste0(prct_LCR_df1, ": ", get(prct_LCR_df1)))
    print(paste0(prct_LCR_df2, ": ", get(prct_LCR_df2)))
    print(paste0(prct_LCR_df3, ": ", get(prct_LCR_df3)))
    
    # get percent 
    prct_not_in_LCR_df1 = paste0(i, "_prct_not_in_LCR_bp_", j, "_df_1")
    prct_not_in_LCR_df2 = paste0(i, "_prct_not_in_LCR_bp_", j, "_df_2")
    prct_not_in_LCR_df3 = paste0(i, "_prct_not_in_LCR_bp_", j, "_df_3")
    
    assign( prct_not_in_LCR_df1, paste0(round(get(n_not_LCR_df1)/get(n_total_df1) * 100, 0), "%"))
    assign( prct_not_in_LCR_df2, paste0(round(get(n_not_LCR_df2)/get(n_total_df2) * 100, 0), "%"))
    assign( prct_not_in_LCR_df3, paste0(round(get(n_not_LCR_df3)/get(n_total_df3) * 100, 0), "%"))
    
    print(paste0(prct_not_in_LCR_df1, ": ", get(prct_not_in_LCR_df1)))
    print(paste0(prct_not_in_LCR_df2, ": ", get(prct_not_in_LCR_df2)))
    print(paste0(prct_not_in_LCR_df3, ": ", get(prct_not_in_LCR_df3)))
  }
}


# Defining lables and sliding window for table S7-S9
Labels_2 = c(rep("rAF_lo"), rep("",2), rep("sAF_hi"), rep("",2), rep("rAF_hi"), rep("",2))
sliding_regions = c(rep(c("Total", "In LCR", "Percentage in LCR"), 3))


rename_columns <- function(data) {
  column_rename_mapping <- c("bp.bp10" = "10bps",
                             "bp.bp20" = "20bps",
                             "bp.bp30" = "30bps",
                             "bp.bp40" = "40bps",
                             "Labels_2" = "Labels",
                             "sliding_regions" = "Sliding Regions")
  
  colnames(data) <- column_rename_mapping[na.omit(colnames(data))]
  return(data)
}




bp10 <- c()
bp20 <- c()
bp30 <- c()
bp40 <- c()
df <- "df_1" #change df depending on which data set you want the table for 

conditions_ordered = c("rAF_lo", "sAF_hi", "rAF_hi")
bp <- list(bp10 = bp10, bp20 = bp20, bp30 = bp30, bp40 = bp40)

information_type <- c('_total_', '_LCR_', '_prct_LCR_')

#Iterates through the base pair ranges, conditions ie sAF_hi / rAF_lo / rAF_hi and if it 
#is the total number The number in LCR and percent in LCR
#a string is made using the iterated information, that string is the variable name from the supplementary table output above

for (j in bp_range) {
  for (i in conditions_ordered) {
    for (z in information_type){
      
      #sets variable_name equal to the variable represented by the string created through the for loops
      variable_name <- get(paste0(i, z, "bp_", j, "_", df))
      
      #sets the bp vector in the bp list equal to itself plus the new variable_name 
      bp[[paste0("bp", j)]] <- c(bp[[paste0("bp", j)]], variable_name)
    }
  }
}



table_S7 = data.frame(Labels_2, sliding_regions, bp$bp10, bp$bp20, bp$bp30, bp$bp40) 


modified_table_S7 <- rename_columns(table_S7)
modified_table_S7
modified_table_S7_df = as.data.frame(modified_table_S7)
setwd("/Users/User_1/Desktop/rAF_project_2/")
fwrite(modified_table_S7_df, "Table S7. The proportion of genomic regions in the gnomAD dataset overlapping with Low Complexity Regions (LCR).csv")



bp10 <- c()
bp20 <- c()
bp30 <- c()
bp40 <- c()
df <- "df_2" #change df depending on which data set you want the table for 

conditions_ordered = c("rAF_lo", "sAF_hi", "rAF_hi")
bp <- list(bp10 = bp10, bp20 = bp20, bp30 = bp30, bp40 = bp40)

information_type <- c('_total_', '_LCR_', '_prct_LCR_')

for (j in bp_range) {
  for (i in conditions_ordered) {
    for (z in information_type){
      variable_name <- get(paste0(i, z, "bp_", j, "_", df))
      
      bp[[paste0("bp", j)]] <- c(bp[[paste0("bp", j)]], variable_name)
    }
  }
}

table_S8 = data.frame(Labels_2, sliding_regions, bp$bp10, bp$bp20, bp$bp30, bp$bp40) 

modified_table_S8 <- rename_columns(table_S8)
modified_table_S8_df = as.data.frame(modified_table_S8)
setwd("/Users/User_1/Desktop/rAF_project_2/")
fwrite(modified_table_S8_df, "Table S8. The proportion of genomic regions in the IGM dataset overlapping with Low Complexity Regions (LCR).csv")

bp10 <- c()
bp20 <- c()
bp30 <- c()
bp40 <- c()
df <- "df_3"

conditions_ordered = c("rAF_lo", "sAF_hi", "rAF_hi")
bp <- list(bp10 = bp10, bp20 = bp20, bp30 = bp30, bp40 = bp40)

information_type <- c('_total_', '_LCR_', '_prct_LCR_')

for (j in bp_range) {
  for (i in conditions_ordered) {
    for (z in information_type){
      variable_name <- get(paste0(i, z, "bp_", j, "_", df))
      
      bp[[paste0("bp", j)]] <- c(bp[[paste0("bp", j)]], variable_name)
    }
  }
}

table_S9 = data.frame(Labels_2, sliding_regions, bp$bp10, bp$bp20, bp$bp30, bp$bp40) 



modified_table_S9 <- rename_columns(table_S9)
modified_table_S9_df = as.data.frame(modified_table_S9)
setwd("/Users/User_1/Desktop/rAF_project_2/")
fwrite(modified_table_S9_df, "Table S9. The proportion of genomic regions in the UK.BB dataset overlapping with Low Complexity Regions (LCR).csv")



for (i in conditions) {
  
  regions = c(rep(paste0(bp_range, " bps"), 2))
  region = c(rep("in LCR", 4), rep("outside of LCR", 4), rep("in LCR", 4), rep("outside of LCR", 4), rep("in LCR", 4), rep("outside of LCR", 4))
  df_labels = c(rep(df_name_1, 8), rep(df_name_2, 8), rep(df_name_3, 8))
  num_indels = c()
  prct_indels = c()
  prct_indels_LCR_only = c()
  
  
  dfs = c("1", "2", "3")
  reg = c("LCR", "not_in_LCR")
  for (j in dfs) {
    for (k in reg) {
      for (l in bp_range) {
        num_var = paste0(i, "_", k, "_bp_", l, "_df_", j)
        num_indels = append(num_indels, get(num_var))
        #print(num_indels)
        
        prct_var = paste0(i, "_prct_", k, "_bp_", l, "_df_", j)
        prct_indels = append(prct_indels, get(prct_var))
        print(prct_indels)
        #this makes the vector that stores only the percent overlap values 
        if (k == "LCR"){
          prct_var_LCR_only = paste0(i, "_prct_", k, "_bp_", l, "_df_", j)
          prct_indels_LCR_only = append(prct_indels_LCR_only, get(prct_var_LCR_only))
          #checks if the length of percent_indels_LCR_only is devisable by 4, if it is, four empty strings are added to the vector, 
          #this vecotr is made to only show the percents attributed to the top bar in the stacked bar chart
          if(length(prct_indels_LCR_only) %% 4 == 0){
            prct_indels_LCR_only = append(prct_indels_LCR_only, c(" ", " ", " ", " "))
          }
        }
      }
    }
  }
  
  graph_df = data.frame(regions, region, num_indels, prct_indels_LCR_only, df_labels)
  
  
  graph = paste0("graph_", i)
  assign(graph, ggplot(graph_df, aes(x = df_labels, y = num_indels, fill = interaction(region, df_labels))) + scale_y_continuous (labels = formatter10000) + geom_bar(stat = "identity", position = "stack") + facet_grid(~regions, switch = "both") + labs(title = paste0(i , " Regions"), y = "Number of Regions (x10000)", x = "regions", fill = "Location")  + scale_fill_grey(labels = c(paste0("Overlap"), paste0 ("Non-Overlap"), paste0("Overlap"), paste0 ("Non-Overlap"),paste0("Overlap"), paste0 ("Non-Overlap"))) + theme(plot.title = element_text(hjust = 0.5, size = 32, margin = margin( 1, 1, 1, 1, "cm")), legend.position = "top", axis.text=element_text(size=26), axis.title=element_text(size=26), legend.title = element_text(size = 26), legend.text = element_text(size = 26), strip.text.x = element_text(size = 26), plot.margin = margin( 1, 1, 1, 1, "cm"), axis.title.x = element_blank(), axis.title.y = element_text(margin = margin( 1, 1, 1, 1, "cm")))) 
  
  if (i == "rAF_lo") {
    assign(graph, get(graph) + labs(title = "rAF-lo Regions") +theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous (labels = formatter10000, limits = c(0,800000)) + geom_bar_text(aes(label = prct_indels_LCR_only), position = "stack", reflow = TRUE, contrast = FALSE, size = 32, padding.y = grid::unit(6, "mm") ,outside = TRUE))
  }
  
  if (i == "sAF_hi"){
    assign(graph, get(graph) + labs(title = "sAF-hi Regions")+theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous (labels = formatter10000, limits = c(0,70000)) + geom_bar_text(aes(label = prct_indels_LCR_only), position = "stack", reflow = TRUE, contrast = FALSE, size = 32, padding.y = grid::unit(-10, "mm") ,outside = TRUE))
  }
  if(i == "rAF_hi"){
    assign(graph, get(graph)+ labs(title = "rAF-hi Regions") + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous (labels = formatter10000, limits = c(0,70000)) + geom_bar_text(aes(label = prct_indels_LCR_only), position = "stack", reflow = TRUE, contrast = FALSE, size = 32, padding.y = grid::unit(-10, "mm") ,outside = TRUE))
    
  }
  
  
}


null_df = data.frame() 
null_figure = ggplot(null_df) + labs(title = "LCR Overlap Analysis") + theme(panel.background = element_rect(fill = 'white', colour = 'white')) + theme (plot.title = element_text(hjust = 0.5, size = 32, margin = margin( 1, 1, 1, 1, "cm")))

graph_summary = ggarrange(print(null_figure), print(graph_rAF_lo), print(graph_sAF_hi), print(graph_rAF_hi), labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, font.label = list(size = 32, color = "black")) + theme(plot.margin = unit(c(2,2,2,2), "cm")) #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

setwd("/Users/User_1/Desktop/rAF_project_2/")
ggsave("Fig4.jpg", width = 80, height = 60, units = c("cm"), dpi = 300)

db = c(rep("gnomad",3), rep("IGM", 3), rep("UK.BB", 3))
type = c(rep(c("rAF_hi", "rAF_lo", "sAF_hi"), 1))
regions_in_LCR = c(rAF_hi_LCR_bp_10_df_1, rAF_lo_LCR_bp_10_df_1, sAF_hi_LCR_bp_10_df_1, rAF_hi_LCR_bp_10_df_2, rAF_lo_LCR_bp_10_df_2, sAF_hi_LCR_bp_10_df_2,rAF_hi_LCR_bp_10_df_3, rAF_lo_LCR_bp_10_df_3, sAF_hi_LCR_bp_10_df_3)
regions_outside_LCR = c(rAF_hi_not_in_LCR_bp_10_df_1, rAF_lo_not_in_LCR_bp_10_df_1, sAF_hi_not_in_LCR_bp_10_df_1, rAF_hi_not_in_LCR_bp_10_df_2, rAF_lo_not_in_LCR_bp_10_df_2, sAF_hi_not_in_LCR_bp_10_df_2, rAF_hi_not_in_LCR_bp_10_df_3, rAF_lo_not_in_LCR_bp_10_df_3, sAF_hi_not_in_LCR_bp_10_df_3)


df = data.frame(db, type, regions_in_LCR, regions_outside_LCR)

fwrite(df, "2023-05-03_p_value_table_LCR.csv")