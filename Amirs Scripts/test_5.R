# Figure 5

setwd("/Users/User_1/Desktop/rAF_project_2")

# Load in packages
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggtext)

#Clear environment: 
rm(list=ls())

# Set your datasets below and read in the files containing all indels ≤ 50bp for each dataset:

## define dfs 
df_name_1 = "gnomAD"
df_name_2 = "IGM"
df_name_3 = "UK.BB"

## Read in dataframe (total)
df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_lt50bp.csv")
df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_lt50bp.csv")
df_3 = fread("/Users/User_1/Desktop/rAF_project_2/UK.BB.exomes.430k.sites_indelsonly_rAF_lt50bp.csv")

#Read in the rAF_hi, rAF_lo and sAF_hi files. 
#rAF_hi indels are sAF ≤ 10^-4 and rAF > 10^-4. 
# rAF_lo indels are sAF ≤ 10^-4 and rAF ≤ 10^-4. 
# sAF_hi indels are sAF > 10^-4.

## Read in rAF_hi 

rAF_hi_bp_10_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_bp_20_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_bp_30_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_bp_40_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

rAF_hi_bp_10_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_bp_20_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_bp_30_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_bp_40_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

rAF_hi_bp_10_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.csv")
rAF_hi_bp_20_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_rAF_hiIndels.lt50bp.csv")
rAF_hi_bp_30_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_rAF_hiIndels.lt50bp.csv")
rAF_hi_bp_40_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_rAF_hiIndels.lt50bp.csv")

## Read in rAF_lo 



rAF_lo_bp_10_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_rAF_loIndels.lt50bp.csv")
rAF_lo_bp_20_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_rAF_loIndels.lt50bp.csv")
rAF_lo_bp_30_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_rAF_loIndels.lt50bp.csv")
rAF_lo_bp_40_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_rAF_loIndels.lt50bp.csv")

rAF_lo_bp_10_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_rAF_loIndels.lt50bp.csv")
rAF_lo_bp_20_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp20_rAF_loIndels.lt50bp.csv")
rAF_lo_bp_30_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp30_rAF_loIndels.lt50bp.csv")
rAF_lo_bp_40_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp40_rAF_loIndels.lt50bp.csv")

rAF_lo_bp_10_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_rAF_loIndels.lt50bp.csv")
rAF_lo_bp_20_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_rAF_loIndels.lt50bp.csv")
rAF_lo_bp_30_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_rAF_loIndels.lt50bp.csv")
rAF_lo_bp_40_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_rAF_loIndels.lt50bp.csv")


## Read in sAF_hi


sAF_hi_bp_10_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_sAF_hiIndels.lt50bp.csv")
sAF_hi_bp_20_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp20_sAF_hiIndels.lt50bp.csv")
sAF_hi_bp_30_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp30_sAF_hiIndels.lt50bp.csv")
sAF_hi_bp_40_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp40_sAF_hiIndels.lt50bp.csv")

sAF_hi_bp_10_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_sAF_hiIndels.lt50bp.csv")
sAF_hi_bp_20_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp20_sAF_hiIndels.lt50bp.csv")
sAF_hi_bp_30_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp30_sAF_hiIndels.lt50bp.csv")
sAF_hi_bp_40_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp40_sAF_hiIndels.lt50bp.csv")

sAF_hi_bp_10_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_sAF_hiIndels.lt50bp.csv")
sAF_hi_bp_20_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp20_sAF_hiIndels.lt50bp.csv")
sAF_hi_bp_30_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp30_sAF_hiIndels.lt50bp.csv")
sAF_hi_bp_40_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp40_sAF_hiIndels.lt50bp.csv")


#Define rAF_lo. 
rAF_lo = (1 * 10^-4)

#Read in ClinVar TSV files: 
clinvar_2 = fread("/Users/User_1/Desktop/rAF_project_2/clinvar38_20230730_onlyClinSig.txt",header = FALSE, quote="")
colnames(clinvar_2) = c("V1", "V2", "V3", "V4", "V10")
clinvar = fread("/Users/User_1/Desktop/rAF_project_2/clinvar_20230318_onlyClinSig.txt", header = FALSE, quote="")
colnames(clinvar) = c("V1", "V2", "V3", "V4", "V10")


#Define ClinVar classfications. 

benign = c('Benign', 'Benign/Likely_benign', 'Benign|_drug_response')
LB = c('Likely_benign', 'Likely_benign|_drug_response|_other',
       'Likely_benign|_other')
conflict = c('Conflicting_interpretations_of_pathogenicity',
             'Conflicting_interpretations_of_pathogenicity|_association',
             'Conflicting_interpretations_of_pathogenicity|_drug_response',
             'Conflicting_interpretations_of_pathogenicity|_drug_response|_other',
             'Conflicting_interpretations_of_pathogenicity|_other',
             'Conflicting_interpretations_of_pathogenicity|_risk_factor',
             'Uncertain_significance',
             'Uncertain_significance|_risk_factor')  #conflict and uncertain significance
LP = c('Likely_pathogenic',
       'Likely_pathogenic|_risk_factor')
patho = c('Pathogenic',
          'Pathogenic/Likely_pathogenic',
          'Pathogenic/Likely_pathogenic|_other', 'Pathogenic|_Affects',
          'Pathogenic|_other', 'Pathogenic|_protective',
          'Pathogenic|_risk_factor')
other = c('Affects', '\\N', 'association',
          'confers_sensitivity', 'drug_response', 'not_provided', 'other',
          'protective', 'risk_factor')

#Define bp ranges for sliding regions. 
bp_range = c("10", "20", "30", "40")

#Give ClinVar file a VarID column. 
clinvar$VarID = paste0(clinvar$V1, "-", clinvar$V2, "-", clinvar$V3, "-", clinvar$V4)
clinvar_2$VarID = paste0(clinvar_2$V1, "-", clinvar_2$V2, "-", clinvar_2$V3, "-", clinvar_2$V4)

# Merge sAF and rAF with ClinVar and keep only the rows that are in ClinVar. Remove any duplicates, if any. 
## Merge sAF and rAF with clinvar and ONLY KEEP THE ROWS THAT ARE IN CLINVAR 

df_1_clinvar = merge(df_1, clinvar, by = "VarID") 
df_2_clinvar = merge(df_2, clinvar, by = "VarID") 
df_3_clinvar = merge(df_3, clinvar_2, by = "VarID") 

## Remove any duplicates, if any 
df_1_clinvar = distinct(df_1_clinvar) 
df_2_clinvar = distinct(df_2_clinvar) 
df_3_clinvar = distinct(df_3_clinvar) 

# Identify rAF_lo sAF indels that are benign/likely benign or pathogenic/likely pathogenic in both datasets. 

## gnomAD 

df_1_B_LB = df_1_clinvar %>% filter(V10 %in% benign | V10 %in% LB) %>% filter (sAF <= rAF_lo) 
df_1_P_LP = df_1_clinvar %>% filter(V10 %in% patho | V10 %in% LP) %>% filter (sAF <= rAF_lo) 


## IGM
df_2_B_LB = df_2_clinvar %>% filter(V10 %in% benign | V10 %in% LB) %>% filter (sAF <= rAF_lo) 
df_2_P_LP = df_2_clinvar %>% filter(V10 %in% patho | V10 %in% LP) %>% filter (sAF <= rAF_lo) 

## UK.BB

df_3_B_LB = df_3_clinvar %>% filter(V10 %in% benign | V10 %in% LB) %>% filter (sAF <= rAF_lo) 
df_3_P_LP = df_3_clinvar %>% filter(V10 %in% patho | V10 %in% LP) %>% filter (sAF <= rAF_lo) 

# Label indels with a log_key. A benign/likely benign indel has a log_key of 0. A pathogenic/likely pathogenic indel has a log key of 1. 

## Give benign/LB a 0 log_key
## Give patho/LP a 1 log_key 
df_1_B_LB$log_key = 0 
df_1_P_LP$log_key = 1 

# #IGM
df_2_B_LB$log_key = 0 
df_2_P_LP$log_key = 1 

#UK.BB
df_3_B_LB$log_key = 0 
df_3_P_LP$log_key = 1 

# Create output dataframes for bargraphs and logistic regression graphs. For the bar graph, we want to zoom into rAF < 0.05 for the 10bp region. 
df_1_output = rbind(df_1_B_LB, df_1_P_LP) %>% select (sAF, rAF_bp10, log_key)

df_2_output = rbind(df_2_B_LB, df_2_P_LP)  %>% select (sAF, rAF_bp10, log_key)

df_3_output = rbind(df_3_B_LB, df_3_P_LP)  %>% select (sAF, rAF_bp10, log_key)

## Mean Square Error 

## Make logistic regression of the rAF/sAF dataframes 
model_rAF_df1  =  glm(formula = log_key ~ rAF_bp10, family = "binomial", data = df_1_output)
model_rAF_df2  = glm(formula = log_key ~ rAF_bp10 , family = "binomial", data = df_2_output)
model_rAF_df3  = glm(formula = log_key ~ rAF_bp10 , family = "binomial", data = df_3_output)

model_sAF_df1  = glm(formula = log_key ~ sAF, family = "binomial", data = df_1_output)
model_sAF_df2  =  glm(formula = log_key ~ sAF, family = "binomial", data = df_2_output)
model_sAF_df3  =  glm(formula = log_key ~ sAF, family = "binomial", data = df_3_output)

## Make dataframe with predicted and actual values 
pred_actual_rAF_df1 = data.frame(pred = predict(model_rAF_df1), actual = df_1_output$rAF_bp10)
pred_actual_rAF_df2 = data.frame(pred = predict(model_rAF_df2), actual = df_2_output$rAF_bp10)
pred_actual_rAF_df3 = data.frame(pred = predict(model_rAF_df3), actual = df_3_output$rAF_bp10)

pred_actual_sAF_df1 = data.frame(pred = predict(model_sAF_df1), actual = df_1_output$sAF)
pred_actual_sAF_df2 = data.frame(pred = predict(model_sAF_df2), actual = df_2_output$sAF)
pred_actual_sAF_df3 = data.frame(pred = predict(model_sAF_df3), actual = df_3_output$sAF)


## Mean square error calculation 
MSE_rAF_df1 = mean((pred_actual_rAF_df1$actual - pred_actual_rAF_df1$pred)^2)
MSE_rAF_df2 = mean((pred_actual_rAF_df2$actual - pred_actual_rAF_df2$pred)^2)
MSE_rAF_df3 = mean((pred_actual_rAF_df3$actual - pred_actual_rAF_df3$pred)^2)
MSE_sAF_df1 = mean((pred_actual_sAF_df1$actual - pred_actual_sAF_df1$pred)^2)
MSE_sAF_df2 = mean((pred_actual_sAF_df2$actual - pred_actual_sAF_df2$pred)^2)
MSE_sAF_df3 = mean((pred_actual_sAF_df3$actual - pred_actual_sAF_df3$pred)^2)

MSE_rAF_df1
MSE_rAF_df2
MSE_rAF_df3
MSE_sAF_df1
MSE_sAF_df2
MSE_sAF_df3


####### For the bargraph, we want to zoom into rAF < 0.05. 
bar_df_1 = df_1_output %>% filter(rAF_bp10 < 0.05)
bar_df_2 = df_2_output %>% filter(rAF_bp10 < 0.05)
bar_df_3 = df_3_output %>% filter(rAF_bp10 < 0.05)

#Create bargraphs.

bargraph_df_1 = ggplot (bar_df_1, aes(x = rAF_bp10 , fill = as.factor(log_key))) + geom_histogram(binwidth = 0.001, position = "dodge") + scale_y_continuous (trans = scales::pseudo_log_trans(), breaks = c(1,100, 1000, 10000)) + scale_fill_manual(labels = c("B/LB", "P/LP"), values = c("#000000", "#989898")) + labs( y = "gnomAD Indels (log)", x = "rAF (10bp Sliding region)", fill = "Classification") +  theme(axis.text=element_text(size=20),  axis.title=element_text(size=20), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.position = c(0.85,0.8), legend.background = element_rect(fill = "white", color = "grey"), panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "white"), plot.margin = margin( 1, 1, 1, 1, "cm"), plot.title = element_text(hjust = 0.5, size = 32, margin = margin( 1, 1, 1, 1, "cm")), axis.title.y = element_text(margin = margin( 0.5, 0.5, 0.5, 0.5, "cm"), angle = 90), axis.title.x = element_text(margin = margin( 1, 1, 1, 1, "cm"))) 


bargraph_df_2 = ggplot (bar_df_2, aes(x = rAF_bp10 , fill = as.factor(log_key))) + geom_histogram(binwidth = 0.001, position = "dodge") + scale_y_continuous (trans = scales::pseudo_log_trans(), breaks = c(1,100, 1000, 10000)) + scale_fill_manual(labels = c("B/LB", "P/LP"), values = c("#000000", "#989898")) + labs( y = "IGM Indels (log)", x = "rAF (10bp Sliding region)", fill = "Classification") +  theme(axis.text=element_text(size=20),  axis.title=element_text(size=20), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.position = c(0.85, 0.8), legend.background = element_rect(fill = "white", color = "grey") , panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "white"), plot.margin = margin( 1, 1, 1, 1, "cm"), plot.title = element_text(hjust = 0.5, size = 32, margin = margin( 1, 1, 1, 1, "cm")), axis.title.y = element_text(margin = margin( 0.5, 0.5, 0.5, 0.5, "cm"), angle = 90), axis.title.x = element_text(margin = margin( 1, 1, 1, 1, "cm")))


bargraph_df_3 = ggplot (bar_df_3, aes(x = rAF_bp10 , fill = as.factor(log_key))) + geom_histogram(binwidth = 0.001, position = "dodge") + scale_y_continuous (trans = scales::pseudo_log_trans(), breaks = c(1,100, 1000, 10000)) + scale_fill_manual(labels = c("B/LB", "P/LP"), values = c("#000000", "#989898")) + labs( y = "UK.BB Indels (log)", x = "rAF (10bp Sliding region)", fill = "Classification") +  theme(axis.text=element_text(size=20),  axis.title=element_text(size=20), legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.position = c(0.85, 0.8), legend.background = element_rect(fill = "white", color = "grey") , panel.background = element_rect(fill = "white", color = "grey"), panel.grid.major = element_line(color = "white"), plot.margin = margin( 1, 1, 1, 1, "cm"), plot.title = element_text(hjust = 0.5, size = 32, margin = margin( 1, 1, 1, 1, "cm")), axis.title.y = element_text(margin = margin( 0.5, 0.5, 0.5, 0.5, "cm"), angle = 90), axis.title.x = element_text(margin = margin( 1, 1, 1, 1, "cm")))


#Create logistic regression graphs with rAF on the x-axis. 

logreg_df_1 = ggplot(df_1_output, aes(x=rAF_bp10, y=log_key)) + geom_point() +
  stat_smooth(method="glm", color="dark grey", se=TRUE,
              method.args = list(family=binomial)) + xlim (0,1.5) + labs( y = "Log Odds (P/LP)", x = "gnomAD rAF") +  theme(axis.text=element_text(size=20),  axis.title=element_text(size=20), plot.margin = margin( 1, 1, 1, 1, "cm"), axis.title.y = element_text(margin = margin( 1, 1, 1, 1, "cm"), angle = 90), axis.title.x = element_text(margin = margin( 1, 1, 1, 1, "cm"))) + geom_text(x = 1.5 ,y = 0.9,inherit.aes = FALSE,label = paste0(""),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + coord_cartesian(clip = "off") + geom_text(x = -0.23 ,y = 1,inherit.aes = FALSE,label = paste0("P/LP" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + geom_text(x = -0.23 ,y = 0,inherit.aes = FALSE,label = paste0("B/LB" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10)


logreg_df_2 = ggplot(df_2_output, aes(x=rAF_bp10, y=log_key)) + geom_point() +
  stat_smooth(method="glm", color="dark grey", se=TRUE,
              method.args = list(family=binomial)) + labs( y = "Log Odds (P/LP)", x = "IGM rAF") +  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), plot.margin = margin( 1, 1, 1, 1, "cm"), axis.title.y = element_text(margin = margin(1, 1, 1, 1, "cm"), angle = 90), axis.title.x = element_text(margin = margin( 1, 1, 1, 1, "cm"))) + geom_text(x = 1.5 ,y = 0.9,inherit.aes = FALSE,label = paste0(""),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + coord_cartesian(clip = "off") + geom_text(x = -0.23 ,y = 1,inherit.aes = FALSE,label = paste0("P/LP" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + geom_text(x = -0.23 ,y = 0,inherit.aes = FALSE,label = paste0("B/LB" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10)


logreg_df_3 = ggplot(df_3_output, aes(x=rAF_bp10, y=log_key)) + geom_point() +
  stat_smooth(method="glm", color="dark grey", se=TRUE,
              method.args = list(family=binomial))+ xlim (0,1.5) + labs( y = "Log Odds (P/LP)", x = "UK.BB rAF") +  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), plot.margin = margin( 1, 1, 1, 1, "cm"), axis.title.y = element_text(margin = margin(1, 1, 1, 1, "cm"), angle = 90), axis.title.x = element_text(margin = margin( 1, 1, 1, 1, "cm"))) + geom_text(x = 1.5 ,y = 0.9,inherit.aes = FALSE,label = paste0(""),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + coord_cartesian(clip = "off") + geom_text(x = -0.23 ,y = 1,inherit.aes = FALSE,label = paste0("P/LP" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + geom_text(x = -0.23 ,y = 0,inherit.aes = FALSE,label = paste0("B/LB" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10)


#Find the logistic regression for sAF 

rAF_logreg_df_1 = ggplot(df_1_output, aes(x=sAF, y=log_key)) + geom_point() +
  stat_smooth(method="glm", color="dark grey", se=TRUE,
              method.args = list(family=binomial)) + xlim (0, 0.0001) + labs( y = "Log Odds (P/LP)", x = "gnomAD sAF") +  theme(axis.text=element_text(size=20),  axis.title=element_text(size=20) , plot.margin = margin( 1, 1, 1, 1, "cm"), axis.title.y = element_text(margin = margin( 1, 1, 1, 1, "cm"), angle = 90), axis.title.x = element_text(margin = margin( 1, 1, 1, 1, "cm"))) + geom_text(x = 0.0001 ,y = 0.9,inherit.aes = FALSE,label = paste0(""),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + coord_cartesian(clip = "off") + geom_text(x = -0.000015 ,y = 1,inherit.aes = FALSE,label = paste0("P/LP" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + geom_text(x =  -0.000015 ,y = 0,inherit.aes = FALSE,label = paste0("B/LB" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10)

rAF_logreg_df_2 = ggplot(df_2_output, aes(x=sAF, y=log_key)) + geom_point() +
  stat_smooth(method="glm", color="dark grey", se=TRUE,
              method.args = list(family=binomial)) + xlim (0, 0.0001) + labs( y = "Log Odds (P/LP)", x = "IGM sAF") +  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), plot.margin = margin( 1, 1, 1, 1, "cm"), axis.title.y = element_text(margin = margin( 1, 1, 1, 1, "cm"), angle = 90), axis.title.x = element_text(margin = margin( 1, 1, 1, 1, "cm"))) + geom_text(x = 0.0001 ,y = 0.9,inherit.aes = FALSE,label = paste0(""),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + coord_cartesian(clip = "off") + geom_text(x =  -0.000015 ,y = 1,inherit.aes = FALSE,label = paste0("P/LP" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + geom_text(x =  -0.000015 ,y = 0,inherit.aes = FALSE,label = paste0("B/LB" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10)


rAF_logreg_df_3 = ggplot(df_3_output, aes(x=sAF, y=log_key)) + geom_point() +
  stat_smooth(method="glm", color="dark grey", se=TRUE,
              method.args = list(family=binomial)) + xlim (0, 0.0001) + labs( y = "Log Odds (P/LP)", x = "UK.BB sAF") +  theme(axis.text=element_text(size=20), axis.title=element_text(size=20), plot.margin = margin( 1, 1, 1, 1, "cm"), axis.title.y = element_text(margin = margin( 1, 1, 1, 1, "cm"), angle = 90), axis.title.x = element_text(margin = margin( 1, 1, 1, 1, "cm"))) + geom_text(x = 0.0001 ,y = 0.9,inherit.aes = FALSE,label = paste0("") ,check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + coord_cartesian(clip = "off") + geom_text(x =  -0.000015 ,y = 1,inherit.aes = FALSE,label = paste0("P/LP" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10) + geom_text(x =  -0.000015 ,y = 0,inherit.aes = FALSE,label = paste0("B/LB" ),check_overlap = TRUE,hjust = 1,fontface = 'bold',size = 10)



graph_summary = ggarrange(print(bargraph_df_1), 
                          print(logreg_df_1), 
                          print(rAF_logreg_df_1), 
                          print(bargraph_df_2), 
                          print(logreg_df_2), 
                          print(rAF_logreg_df_2), 
                          print(bargraph_df_3), 
                          print(logreg_df_3), 
                          print(rAF_logreg_df_3), 
                          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
                          ncol = 3, nrow = 3, font.label = list(size = 32, color = "black")) + theme(plot.margin = unit(c(2,2,2,2), "cm"))


graph_summary 

ggsave("Fig5_NEW.jpg", width = 80, height = 60, units = c("cm"), dpi = 300)


conditions = c("rAF_hi", "rAF_lo", "sAF_hi")

categories = c("benign", "LB", "VUS/other", "LP", "patho" )

for (i in conditions ) {
  for (j in categories) {
    ## get temp variable for the condition 
    condition_df1 = paste0(i, "_bp_10_df_1")
    condition_df2 = paste0(i, "_bp_10_df_2")
    condition_df3 = paste0(i, "_bp_10_df_3")
    
    ## set a temp variable for a clinvar merge 
    clinvar_cond_df1 = paste0("clinvar_", i, "_bp_10_df_1")
    clinvar_cond_df2 = paste0("clinvar_", i, "_bp_10_df_2")
    clinvar_cond_df3 = paste0("clinvar_", i, "_bp_10_df_3")
    
    ## merge with clinvar
    assign (clinvar_cond_df1, merge(get(condition_df1), clinvar, by = "VarID"))
    assign (clinvar_cond_df2, merge(get(condition_df2), clinvar, by = "VarID"))
    ## we think that there is smth going on with the clinvar merging here. 
    assign (clinvar_cond_df3, merge(get(condition_df3), clinvar_2, by = "VarID"))
    
    ## get a dataframe only for the category (ie. benign, LB, conflict, LP, patho, other)
    if (j != "VUS/other") {
      ##  set a temp variable for a clinvar merge
      category_cond_df1 = paste0(j , "_", i, "_bp_10_df1")
      category_cond_df2 = paste0(j , "_", i, "_bp_10_df2")
      category_cond_df3 = paste0(j , "_", i, "_bp_10_df3")
      
      #assign temp variable
      assign (category_cond_df1, get(clinvar_cond_df1) %>% filter(V10 %in% get(j)))
      assign (category_cond_df2, get(clinvar_cond_df2) %>% filter(V10 %in% get(j)))
      assign (category_cond_df3, get(clinvar_cond_df3) %>% filter(V10 %in% get(j)))
    }
    
    else {
      ##  set a temp variable for a clinvar merge
      category_cond_df1 = paste0("VUS_other_", i, "_bp_10_df1")
      category_cond_df2 = paste0("VUS_other_", i, "_bp_10_df2")
      category_cond_df3 = paste0("VUS_other_", i, "_bp_10_df3")
      
      #assign temp variable
      assign (category_cond_df1, get(clinvar_cond_df1) %>% filter(V10 %in% other | V10 %in% conflict))
      assign (category_cond_df2, get(clinvar_cond_df2) %>% filter(V10 %in% other | V10 %in% conflict))
      assign (category_cond_df3, get(clinvar_cond_df3) %>% filter(V10 %in% other | V10 %in% conflict))
    }
    
  }
}


## make the df for each category for db 1 

sAF_hi_categories_df1 = c(nrow(benign_sAF_hi_bp_10_df1), nrow(LB_sAF_hi_bp_10_df1), nrow(VUS_other_sAF_hi_bp_10_df1), nrow(LP_sAF_hi_bp_10_df1), nrow(patho_sAF_hi_bp_10_df1))

rAF_hi_categories_df1 = c(nrow(benign_rAF_hi_bp_10_df1), nrow(LB_rAF_hi_bp_10_df1), nrow(VUS_other_rAF_hi_bp_10_df1), nrow(LP_rAF_hi_bp_10_df1), nrow(patho_rAF_hi_bp_10_df1))

rAF_lo_categories_df1 = c(nrow(benign_rAF_lo_bp_10_df1), nrow(LB_rAF_lo_bp_10_df1), nrow(VUS_other_rAF_lo_bp_10_df1), nrow(LP_rAF_lo_bp_10_df1), nrow(patho_rAF_lo_bp_10_df1))

categories_df1 = data.frame(categories, sAF_hi_categories_df1, rAF_hi_categories_df1, rAF_lo_categories_df1)


## make the df for each category for db 2
sAF_hi_categories_df2 = c(nrow(benign_sAF_hi_bp_10_df2), nrow(LB_sAF_hi_bp_10_df2), nrow(VUS_other_sAF_hi_bp_10_df2), nrow(LP_sAF_hi_bp_10_df2), nrow(patho_sAF_hi_bp_10_df2))

rAF_hi_categories_df2 = c(nrow(benign_rAF_hi_bp_10_df2), nrow(LB_rAF_hi_bp_10_df2), nrow(VUS_other_rAF_hi_bp_10_df2), nrow(LP_rAF_hi_bp_10_df2), nrow(patho_rAF_hi_bp_10_df2))

rAF_lo_categories_df2 = c(nrow(benign_rAF_lo_bp_10_df2), nrow(LB_rAF_lo_bp_10_df2), nrow(VUS_other_rAF_lo_bp_10_df2), nrow(LP_rAF_lo_bp_10_df2), nrow(patho_rAF_lo_bp_10_df2))

categories_df2 = data.frame(categories, sAF_hi_categories_df2, rAF_hi_categories_df2, rAF_lo_categories_df2)

## make the df for each category for db 3
sAF_hi_categories_df3 = c(nrow(benign_sAF_hi_bp_10_df3), nrow(LB_sAF_hi_bp_10_df3), nrow(VUS_other_sAF_hi_bp_10_df3), nrow(LP_sAF_hi_bp_10_df3), nrow(patho_sAF_hi_bp_10_df3))

rAF_hi_categories_df3 = c(nrow(benign_rAF_hi_bp_10_df3), nrow(LB_rAF_hi_bp_10_df3), nrow(VUS_other_rAF_hi_bp_10_df3), nrow(LP_rAF_hi_bp_10_df3), nrow(patho_rAF_hi_bp_10_df3))

rAF_lo_categories_df3 = c(nrow(benign_rAF_lo_bp_10_df3), nrow(LB_rAF_lo_bp_10_df3), nrow(VUS_other_rAF_lo_bp_10_df3), nrow(LP_rAF_lo_bp_10_df3), nrow(patho_rAF_lo_bp_10_df3))

categories_df3 = data.frame(categories, sAF_hi_categories_df3, rAF_hi_categories_df3, rAF_lo_categories_df3)

categories_df1
categories_df2
categories_df3

fwrite(categories_df1, "2023-05-08_pvalues_clinvar_gnomad.csv")
fwrite(categories_df2, "2023-05-08_pvalues_clinvar_igm.csv")
fwrite(categories_df3, "2023-05-08_pvalues_clinvar_UK.BB.csv") #may have to change name scheme later 


# Make Table S9/S10

#Find the percentage of rAF_hi indels for each bp region for each of the following categories: Benign/likely benign (LB), pathogenic (patho) and likey pathogenic (LP) 

cat = c("benign", "LB", "patho", "LP")

## Find rAF_lo 
for (i in cat) {
  ## set var 
  cat_1 = paste0(i, "_rAF_lo_sAF_df1")
  cat_2 = paste0(i, "_rAF_lo_sAF_df2")
  cat_3 = paste0(i, "_rAF_lo_sAF_df3")
  
  #assign
  assign(cat_1, df_1_clinvar %>% filter(sAF <= rAF_lo) %>% filter(V10 %in% get(i)))
  assign(cat_2, df_2_clinvar %>% filter(sAF <= rAF_lo) %>% filter(V10 %in% get(i)))
  assign(cat_3, df_3_clinvar %>% filter(sAF <= rAF_lo) %>% filter(V10 %in% get(i)))
}

## Find percentage
for (i in bp_range) {
  for (j in cat) {
    ## get denominator (rAF_lo)
    rAF_lo_sAF_df1 = paste0(j, "_rAF_lo_sAF_df1")
    rAF_lo_sAF_df2 = paste0(j, "_rAF_lo_sAF_df2")
    rAF_lo_sAF_df3 = paste0(j, "_rAF_lo_sAF_df3")
    
    ## set percent output 
    prct_rAF_hi_1 = paste0("prct_bp_", i, "_", j, "_df1")
    prct_rAF_hi_2 = paste0("prct_bp_", i, "_", j, "_df2")
    prct_rAF_hi_3 = paste0("prct_bp_", i, "_", j, "_df3")
    
    #set total output
    total_rAF_hi_1 = paste0("total_bp_", i, "_", j, "_df1")
    total_rAF_hi_2 = paste0("total_bp_", i, "_", j, "_df2")
    total_rAF_hi_3 = paste0("total_bp_", i, "_", j, "_df3")
    ## set rAF region 
    temp_region = paste0("rAF_bp", i)
    
    ## assign percentage 
    assign (prct_rAF_hi_1, paste0(round(nrow(get(rAF_lo_sAF_df1) %>% filter(get(temp_region) > rAF_lo)) / nrow(get(rAF_lo_sAF_df1)) * 100 , 0), "%"))
    assign(total_rAF_hi_1, paste0(nrow(get(rAF_lo_sAF_df1) %>% filter(get(temp_region) > rAF_lo))))
    
    assign (prct_rAF_hi_2, paste0(round(nrow(get(rAF_lo_sAF_df2) %>% filter(get(temp_region) > rAF_lo)) / nrow(get(rAF_lo_sAF_df2)) * 100 , 0), "%"))
    assign(total_rAF_hi_2, paste0(nrow(get(rAF_lo_sAF_df2) %>% filter(get(temp_region) > rAF_lo))))
    
    assign (prct_rAF_hi_3, paste0(round(nrow(get(rAF_lo_sAF_df3) %>% filter(get(temp_region) > rAF_lo)) / nrow(get(rAF_lo_sAF_df3)) * 100 , 0), "%"))
    assign(total_rAF_hi_3, paste0(nrow(get(rAF_lo_sAF_df3) %>% filter(get(temp_region) > rAF_lo))))
    print(total_rAF_hi_3)
    
  }
}


## make df that show the percentage of rAF_hi indels for each categories for dataframe 1 
rAF_lo_sAF_df1 = c(nrow(benign_rAF_lo_sAF_df1), nrow(LB_rAF_lo_sAF_df1), nrow(patho_rAF_lo_sAF_df1), nrow(LP_rAF_lo_sAF_df1))
prct_rAF_hi_10bp_df1 = c(paste0(total_bp_10_benign_df1," (",prct_bp_10_benign_df1,")"), paste0(total_bp_10_LB_df1," (",prct_bp_10_LB_df1,")"), paste0(total_bp_10_patho_df1," (",prct_bp_10_patho_df1,")"), paste0(total_bp_10_LP_df1," (",prct_bp_10_LP_df1,")"))
prct_rAF_hi_20bp_df1 = c(paste0(total_bp_20_benign_df1," (",prct_bp_20_benign_df1,")"), paste0(total_bp_20_LB_df1," (",prct_bp_20_LB_df1,")"), paste0(total_bp_20_patho_df1," (",prct_bp_20_patho_df1,")"), paste0(total_bp_20_LP_df1," (",prct_bp_20_LP_df1,")"))
prct_rAF_hi_30bp_df1 = c(paste0(total_bp_30_benign_df1," (",prct_bp_30_benign_df1,")"), paste0(total_bp_30_LB_df1," (",prct_bp_30_LB_df1,")"), paste0(total_bp_30_patho_df1," (",prct_bp_30_patho_df1,")"), paste0(total_bp_30_LP_df1," (",prct_bp_30_LP_df1,")"))
prct_rAF_hi_40bp_df1 = c(paste0(total_bp_40_benign_df1," (",prct_bp_40_benign_df1,")"), paste0(total_bp_40_LB_df1," (",prct_bp_40_LB_df1,")"), paste0(total_bp_40_patho_df1," (",prct_bp_40_patho_df1,")"), paste0(total_bp_40_LP_df1," (",prct_bp_40_LP_df1,")"))

prct_rAF_hi_cat_df_1 = data.frame (cat, rAF_lo_sAF_df1, prct_rAF_hi_10bp_df1, prct_rAF_hi_20bp_df1, prct_rAF_hi_30bp_df1, prct_rAF_hi_40bp_df1)


## make df that show the percentage of rAF_hi indels for each categories for dataframe 2 
rAF_lo_sAF_df2 = c(nrow(benign_rAF_lo_sAF_df2), nrow(LB_rAF_lo_sAF_df2), nrow(patho_rAF_lo_bp_10_df2), nrow(LP_rAF_lo_sAF_df2))
prct_rAF_hi_10bp_df2 = c(paste0(total_bp_10_benign_df2," (",prct_bp_10_benign_df2,")"), paste0(total_bp_10_LB_df2," (",prct_bp_10_LB_df2,")"), paste0(total_bp_10_patho_df2," (",prct_bp_10_patho_df2,")"), paste0(total_bp_10_LP_df2," (",prct_bp_10_LP_df2,")"))
prct_rAF_hi_20bp_df2 = c(paste0(total_bp_20_benign_df2," (",prct_bp_20_benign_df2,")"), paste0(total_bp_20_LB_df2," (",prct_bp_20_LB_df2,")"), paste0(total_bp_20_patho_df2," (",prct_bp_20_patho_df2,")"), paste0(total_bp_20_LP_df2," (",prct_bp_20_LP_df2,")"))
prct_rAF_hi_30bp_df2 = c(paste0(total_bp_30_benign_df2," (",prct_bp_30_benign_df2,")"), paste0(total_bp_30_LB_df2," (",prct_bp_30_LB_df2,")"), paste0(total_bp_30_patho_df2," (",prct_bp_30_patho_df2,")"), paste0(total_bp_30_LP_df2," (",prct_bp_30_LP_df2,")"))
prct_rAF_hi_40bp_df2 = c(paste0(total_bp_40_benign_df2," (",prct_bp_40_benign_df2,")"), paste0(total_bp_40_LB_df2," (",prct_bp_40_LB_df2,")"), paste0(total_bp_40_patho_df2," (",prct_bp_40_patho_df2,")"), paste0(total_bp_40_LP_df2," (",prct_bp_40_LP_df2,")"))

prct_rAF_hi_cat_df_2 = data.frame (cat, rAF_lo_sAF_df2, prct_rAF_hi_10bp_df2, prct_rAF_hi_20bp_df2, prct_rAF_hi_30bp_df2, prct_rAF_hi_40bp_df2)

## make df that show the percentage of rAF_hi indels for each categories for dataframe 2 
rAF_lo_sAF_df3 = c(nrow(benign_rAF_lo_sAF_df3), nrow(LB_rAF_lo_sAF_df3), nrow(patho_rAF_lo_bp_10_df3), nrow(LP_rAF_lo_sAF_df3))  
prct_rAF_hi_10bp_df3 = c(paste0(total_bp_10_benign_df3," (",prct_bp_10_benign_df3,")"), paste0(total_bp_10_LB_df3," (",prct_bp_10_LB_df3,")"), paste0(total_bp_10_patho_df3," (",prct_bp_10_patho_df3,")"), paste0(total_bp_10_LP_df3," (",prct_bp_10_LP_df3,")"))
prct_rAF_hi_20bp_df3 = c(paste0(total_bp_20_benign_df3," (",prct_bp_20_benign_df3,")"), paste0(total_bp_20_LB_df3," (",prct_bp_20_LB_df3,")"), paste0(total_bp_20_patho_df3," (",prct_bp_20_patho_df3,")"), paste0(total_bp_20_LP_df3," (",prct_bp_20_LP_df3,")"))
prct_rAF_hi_30bp_df3 = c(paste0(total_bp_30_benign_df3," (",prct_bp_30_benign_df3,")"), paste0(total_bp_30_LB_df3," (",prct_bp_30_LB_df3,")"), paste0(total_bp_30_patho_df3," (",prct_bp_30_patho_df3,")"), paste0(total_bp_30_LP_df3," (",prct_bp_30_LP_df3,")"))
prct_rAF_hi_40bp_df3 = c(paste0(total_bp_40_benign_df3," (",prct_bp_40_benign_df3,")"), paste0(total_bp_40_LB_df3," (",prct_bp_40_LB_df3,")"), paste0(total_bp_40_patho_df3," (",prct_bp_40_patho_df3,")"), paste0(total_bp_40_LP_df3," (",prct_bp_40_LP_df3,")"))

prct_rAF_hi_cat_df_3 = data.frame (cat, rAF_lo_sAF_df3, prct_rAF_hi_10bp_df3, prct_rAF_hi_20bp_df3, prct_rAF_hi_30bp_df3, prct_rAF_hi_40bp_df3)

prct_rAF_hi_cat_df_1
prct_rAF_hi_cat_df_2
prct_rAF_hi_cat_df_3



## Find the number of rAF_hi by sAF P/LP indels that have a rAF > 1% 


for (i in bp_range) {
  
  temp_df1 = paste0("bp_", i,"_P_LP_rAF_gt_1_df1")
  temp_df2 = paste0("bp_", i,"_P_LP_rAF_gt_1_df2")
  temp_df3 = paste0("bp_", i,"_P_LP_rAF_gt_1_df3")
  
  ## set rAF region 
  temp_region = paste0("rAF_bp", i)
  
  assign(temp_df1, nrow(df_1_clinvar %>% filter (sAF <= rAF_lo & get(temp_region) > rAF_lo) %>%  filter(V10 %in% patho | V10 %in% LP) %>% filter(get(temp_region) > 0.01)))
  assign(temp_df2, nrow(df_2_clinvar %>% filter (sAF <= rAF_lo & get(temp_region) > rAF_lo) %>% filter(V10 %in% patho | V10 %in% LP) %>% filter(get(temp_region) > 0.01)))
  assign(temp_df3, nrow(df_3_clinvar %>% filter (sAF <= rAF_lo & get(temp_region) > rAF_lo) %>% filter(V10 %in% patho | V10 %in% LP) %>% filter(get(temp_region) > 0.01)))
}

P_LP_rAF_gt_1_df1 = c(bp_10_P_LP_rAF_gt_1_df1, bp_20_P_LP_rAF_gt_1_df1, bp_30_P_LP_rAF_gt_1_df1, bp_40_P_LP_rAF_gt_1_df1)
P_LP_rAF_gt_1_df2 = c(bp_10_P_LP_rAF_gt_1_df2, bp_20_P_LP_rAF_gt_1_df2, bp_30_P_LP_rAF_gt_1_df2, bp_40_P_LP_rAF_gt_1_df2)
P_LP_rAF_gt_1_df3 = c(bp_10_P_LP_rAF_gt_1_df3, bp_20_P_LP_rAF_gt_1_df3, bp_30_P_LP_rAF_gt_1_df3, bp_40_P_LP_rAF_gt_1_df3)

dataframe_P_LP_gt_1 = data.frame(bp_range, P_LP_rAF_gt_1_df1, P_LP_rAF_gt_1_df2, P_LP_rAF_gt_1_df3)


dataframe_P_LP_gt_1

#fwrite(dataframe_P_LP_gt_1, "New_SupTable_s9.csv")

## Supplementary Excel 1 and Excel 2 


#Make supplementary excel files that show the number of pathogenic/likely pathogenic and benign/likely benign indels there are for each rAF value (rounded to the 6th decimal place) for both datasets. 

### Supplementary excel 1: for database 1 
df_1_B_LB$rAF_bp10 = round(df_1_B_LB$rAF_bp10, 6 )
df_1_P_LP$rAF_bp10 = round(df_1_P_LP$rAF_bp10, 6 )


df_1_freq_B_LB =  df_1_B_LB %>% group_by(rAF_bp10) %>% summarise(count_B_LB=n())
df_1_freq_P_LP =  df_1_P_LP %>% group_by(rAF_bp10) %>% summarise(count_P_LP=n())

df_1_freq = merge(df_1_freq_B_LB, df_1_freq_P_LP, by = "rAF_bp10", all = TRUE)
df_1_freq[is.na(df_1_freq)] = 0

df_1_freq

fwrite(df_1_freq, "supplementary_excel_1.csv")

### Supplementary excel 1: for database 2
df_2_B_LB$rAF_bp10 = round(df_2_B_LB$rAF_bp10, 6 )
df_2_P_LP$rAF_bp10 = round(df_2_P_LP$rAF_bp10, 6 )


df_2_freq_B_LB =  df_2_B_LB %>% group_by(rAF_bp10) %>% summarise(count_B_LB=n())
df_2_freq_P_LP =  df_2_P_LP %>% group_by(rAF_bp10) %>% summarise(count_P_LP=n())

df_2_freq = merge(df_2_freq_B_LB, df_2_freq_P_LP, by = "rAF_bp10", all = TRUE)
df_2_freq[is.na(df_2_freq)] = 0

df_2_freq

fwrite(df_2_freq, "supplementary_excel_2.csv")


### Supplementary excel 1: for database 3
df_3_B_LB$rAF_bp10 = round(df_3_B_LB$rAF_bp10, 6 )
df_3_P_LP$rAF_bp10 = round(df_3_P_LP$rAF_bp10, 6 )


df_3_freq_B_LB =  df_3_B_LB %>% group_by(rAF_bp10) %>% summarise(count_B_LB=n())
df_3_freq_P_LP =  df_3_P_LP %>% group_by(rAF_bp10) %>% summarise(count_P_LP=n())

df_3_freq = merge(df_3_freq_B_LB, df_3_freq_P_LP, by = "rAF_bp10", all = TRUE)
df_3_freq[is.na(df_3_freq)] = 0

df_3_freq

fwrite(df_3_freq, "supplementary_excel_3.csv")

