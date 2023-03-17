## Sandy Yang 
## February 13, 2023 
## Figure 1 

###################################################
############## load packages ######################
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
library(ggplot2)
install.packages("ggbreak")
library(ggbreak)

## set rare threshold 
rare = 10^-4 

## Read in files 


data = fread("2023-03-06_gnomAD_sAF_rAF_lte_50bp.csv")

print("Finished loading in data.")

dataset = "gnomAD"

## Order by 
  ## chr (ascending)
  ## pos (ascending)
data = data[with(data, order(CHR, POS)), ]

## Needs: 
  ## 1) Make a new col to calculate the length by looking at alternate and ref positions (this abs value should be <= 50) to find indels < 50bp in length
  ## 2) If it is an insertion: length of ALT > length of REF
  ## 3) If it is a deletion: length of ALT < length of REF 
data = data %>% mutate(
  type = case_when(nchar(data$ALT) > nchar(data$REF) ~ "insertion", TRUE ~ "deletion"), 
  indel_length = nchar(data$ALT) - nchar(data$REF))

## Things we need to include: 
  ## all indels in a bin with AT LEAST 1 rare sAF/common rAF indel regardless of rare (means NOT ONLY suspicious indels are kept)
  ## you can tell if an indel is in the same bin if they have the same CHR and also same groupID number for that bp 
  ## If a bin has at least one suspicious indel (rare sAF/common rAF), all the indels in this bin are kept regardless of whether they are also rare sAF/common rAF as long as one indel in that bin is because this helps to maintain the integrity of the bin
  ## if we removed these indels, then indels that are in a bin will no longer be close to each other because the position where you begin to count from will change if you disregard these indels 

## To find the bin with at least one suspicious indel:
## filter for suspicious indels
## get the CHR and groupID pairs 
## Get that inofrmation and then use it for filtering big dataset 
bins_with_at_least_one_sus_indel_bp5 = distinct(data %>% filter(data$sAF <= rare & data$rAF5 > rare ) %>% select ("CHR", "GROUPID_bp5")) #9019
bins_with_at_least_one_sus_indel_bp10 = distinct(data %>% filter(data$sAF <= rare & data$rAF10 > rare ) %>% select ("CHR", "GROUPID_bp10")) #10478
bins_with_at_least_one_sus_indel_bp15 = distinct(data %>% filter(data$sAF <= rare & data$rAF15 > rare ) %>% select ("CHR", "GROUPID_bp15")) #11600
bins_with_at_least_one_sus_indel_bp20 = distinct(data %>% filter(data$sAF <= rare & data$rAF20 > rare ) %>% select ("CHR", "GROUPID_bp20")) #12508

## subset by keeping indels <=50bp in length, and only keep bins with more than one indel in it, and only keep bins with at least one suspicious indel 
## subset
  ## keep bins with at least one sus indel 
  ## keep indels <= 50bp in length 
  ## keep bins with more than one indel in it 

#42865
subset_bp5 = data %>% filter(paste0(data$CHR, "-", data$GROUPID_bp5) %in% paste0(bins_with_at_least_one_sus_indel_bp5$CHR, "-", bins_with_at_least_one_sus_indel_bp5$GROUPID_bp5)) %>% filter(abs(indel_length) <= 50) %>% group_by(CHR, GROUPID_bp5) %>% filter(n()>1)

#55955
subset_bp10 = data %>% filter(paste0(data$CHR, "-", data$GROUPID_bp10) %in% paste0(bins_with_at_least_one_sus_indel_bp10$CHR, "-", bins_with_at_least_one_sus_indel_bp10$GROUPID_bp10)) %>% filter(abs(indel_length) <= 50) %>% group_by(CHR, GROUPID_bp10) %>% filter(n()>1)

#66198
subset_bp15 = data %>% filter(paste0(data$CHR, "-", data$GROUPID_bp15) %in% paste0(bins_with_at_least_one_sus_indel_bp15$CHR, "-", bins_with_at_least_one_sus_indel_bp15$GROUPID_bp15)) %>% filter(abs(indel_length) <= 50) %>% group_by(CHR, GROUPID_bp15) %>% filter(n()>1)

#75013
subset_bp20 = data %>% filter(paste0(data$CHR, "-", data$GROUPID_bp20) %in% paste0(bins_with_at_least_one_sus_indel_bp20$CHR, "-", bins_with_at_least_one_sus_indel_bp20$GROUPID_bp20)) %>% filter(abs(indel_length) <= 50) %>% group_by(CHR, GROUPID_bp20) %>% filter(n()>1)



## calculate bin length 
## if ending indel is insertion: position of ending indel - position of starting indel + length of ending indel if ending indel is an insertion
## If ending indel is deletion: position of ending indel - position of starting indel 

## find the start indel in that group 
subset_bp5 = subset_bp5 %>% group_by(CHR, GROUPID_bp5) %>% mutate (start_indel = min (POS)) %>% ungroup() 
subset_bp10 = subset_bp10 %>% group_by(CHR, GROUPID_bp10) %>% mutate (start_indel = min (POS)) %>% ungroup() 
subset_bp15 = subset_bp15 %>% group_by(CHR, GROUPID_bp15) %>% mutate (start_indel = min (POS)) %>% ungroup() 
subset_bp20 = subset_bp20 %>% group_by(CHR, GROUPID_bp20) %>% mutate (start_indel = min (POS)) %>% ungroup() 

## find the end indel in that group 
subset_bp5 = subset_bp5 %>% group_by(CHR, GROUPID_bp5) %>% mutate (end_indel = max (POS)) %>% ungroup() 
subset_bp10 = subset_bp10 %>% group_by(CHR, GROUPID_bp10) %>% mutate (end_indel = max (POS)) %>% ungroup() 
subset_bp15 = subset_bp15 %>% group_by(CHR, GROUPID_bp15) %>% mutate (end_indel = max (POS)) %>% ungroup() 
subset_bp20 = subset_bp20 %>% group_by(CHR, GROUPID_bp20) %>% mutate (end_indel = max (POS)) %>% ungroup() 

## calculate the bin length 
subset_bp5$bin_len_window_10 = subset_bp5$end_indel - subset_bp5$start_indel
subset_bp10$bin_len_window_20 = subset_bp10$end_indel - subset_bp10$start_indel
subset_bp15$bin_len_window_30 = subset_bp15$end_indel - subset_bp15$start_indel
subset_bp20$bin_len_window_40 = subset_bp20$end_indel - subset_bp20$start_indel

## Get a chart with the number of bins with that bin length 
bins_window_10 = distinct(subset_bp5 %>% select ("CHR", "GROUPID_bp5", "bin_len_window_10")) 
bins_window_20 = distinct(subset_bp10 %>% select ("CHR", "GROUPID_bp10", "bin_len_window_20")) 
bins_window_30 = distinct(subset_bp15 %>% select ("CHR", "GROUPID_bp15", "bin_len_window_30")) 
bins_window_40 = distinct(subset_bp20 %>% select ("CHR", "GROUPID_bp20", "bin_len_window_40")) 

## Make graphs
#window 10 
ggplot(bins_window_10, aes(x = bin_len_window_10)) +
  geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
  geom_density(fill="grey", alpha = .5) + 
  scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) + 
  scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 850, 900, 950, 1000)) +
  ylim(0, 0.3) + 
  labs(title = paste0(dataset, ": Density Plots: Bin Lengths in bp Window 10"), x = "Bin Lengths (bps)")


#window 20
ggplot(bins_window_20, aes(x = bin_len_window_20)) +
  geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
  geom_density(fill="grey", alpha = .5) + 
  scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) + 
  scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 850, 900, 950, 1000)) +
  ylim(0, 0.3) + 
  labs(title = paste0(dataset, ": Density Plots: Bin Lengths in bp Window 20"), x = "Bin Lengths (bps)")

#window 30

ggplot(bins_window_30, aes(x = bin_len_window_30)) +
  geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
  geom_density(fill="grey", alpha = .5) + 
  scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) + 
  scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 850, 900, 950, 1000)) +
  ylim(0, 0.3) + 
  labs(title = paste0(dataset, ": Density Plots: Bin Lengths in bp Window 30"), x = "Bin Lengths (bps)")


#window 40
ggplot(bins_window_40, aes(x = bin_len_window_40)) +
  geom_histogram(aes(y = ..density..), colour= "black", fill = "white", binwidth = 1) +
  geom_density(fill="grey", alpha = .5) + 
  scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) + 
  scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 850, 900, 950, 1000)) +
  ylim(0, 0.3) + 
  labs(title = paste0(dataset, ": Density Plots: Bin Lengths in bp Window 40"), x = "Bin Lengths (bps)")


## combine graphs together 
bins_window_10$window <- "10bp_window"
bins_window_20$window <- "20bp_window"
bins_window_30$window <- "30bp_window"
bins_window_40$window <- "40bp_window"

add_window10 = bins_window_10 %>% select ("bin_len_window_10", "window")
add_window20 = bins_window_20 %>% select ("bin_len_window_20", "window")
add_window30 = bins_window_30 %>% select ("bin_len_window_30", "window")
add_window40 = bins_window_40 %>% select ("bin_len_window_40", "window")

colnames(add_window10) = c("bin_len", "window")
colnames(add_window20) = c("bin_len", "window")
colnames(add_window30) = c("bin_len", "window")
colnames(add_window40) = c("bin_len", "window")

combined = rbind (add_window10, add_window20, add_window30, add_window40)

ggplot(combined, aes(x = bin_len, group = window, color= window, fill = window)) +
  geom_histogram(aes(y = ..density.., fill = window), position = "identity", alpha = 0.7, color = "black", binwidth = 1) +
  geom_density(aes(fill= window), alpha = .5) + 
  scale_x_continuous(breaks =seq(from = 0, to = 1000, by = 20), limits = c(0, 1000)) + 
  scale_x_break(breaks = c(80, 800), scales = 0.25, ticklabels = c(800, 850, 900, 950, 1000)) +
  ylim(0, 0.3) + 
  labs(title = paste0(dataset, ": Density Plots"), x = "Bin Lengths (bps)")


  
#####################################
# For written Results section # 
#####################################
  
# What percentage of bins have < 10bp in length in each window? 
count_window10 = bins_window_10 %>% group_by(bin_len_window_10) %>% count(bin_len_window_10)
count_window20 = bins_window_20 %>% group_by(bin_len_window_20) %>% count(bin_len_window_20)
count_window30 = bins_window_30 %>% group_by(bin_len_window_30) %>% count(bin_len_window_30)
count_window40 = bins_window_40 %>% group_by(bin_len_window_40) %>% count(bin_len_window_40)

percent_lt_10bp_window_10 = sum((count_window10 %>% filter(bin_len_window_10 < 10))$n) / sum(count_window10$n) * 100 
percent_lt_20bp_window_20 = sum((count_window20 %>% filter(bin_len_window_20 < 20))$n) / sum(count_window20$n) * 100 
percent_lt_30bp_window_30 = sum((count_window30 %>% filter(bin_len_window_30 < 30))$n) / sum(count_window30$n) * 100 
percent_lt_40bp_window_40 = sum((count_window40 %>% filter(bin_len_window_40 < 40))$n) / sum(count_window40$n) * 100 





