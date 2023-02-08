## Author: Sandy Yang
## December 19, 2022 
## This code calculates the rAF values
## Filter by chr and groupID (iterating over all four types of windows) and then find the mean AN 
## calculate the rAF by finding the AC/meanAN
## for gnomAD data: make sure to multiply by 2 because we are using the threshold 1 * 10^-4

###################################################
############## load packages ######################
library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)
####################################################
##################### set rare #####################
rare = 1 * 10^-4
####################################################
################# set database #####################
database = "gnomAD" #or IGM
####################################################
############## read in the files ###################
#For gnomAD
if (database == "gnomAD") {
  bin_ID = fread("gnomad.exomes.r2.1.1.sites_indelsonly_allele_counts.csv")
  indels_test = fread("20220926_gnomad_srAF_double_rare_threshold_10_4.csv")
  indels = fread("gnomad_indels_with_sAF.csv")
}

#For IGM
if (database == "IGM") {
  bin_ID = fread("20220906_IGM__VarID_AC_AN_no_duplic_sorted_header.csv")
  indels_test = fread("20220909_IGM_srAF_sorted.csv")
  indels = fread("IGM_indels_with_sAF.csv")
}
####################################################
################## set VarID #######################
## VarID
bin_ID$VarID = paste0(bin_ID$CHR, "-", bin_ID$POS, "-", bin_ID$REF, "-", bin_ID$ALT)
indels$VarID = paste0(indels$CHR, "-", indels$POS, "-", indels$REF, "-", indels$ALT)
####################################################
######## create data set with all info #############
bin_ID = bin_ID[ , c("GROUPID_bp5", "GROUPID_bp10", "GROUPID_bp15", "GROUPID_bp20", "AlleleCount_bp5", "AlleleCount_bp10", "AlleleCount_bp15", "AlleleCount_bp20", "VarID")]
## left_join 
bin_indels = left_join(indels, bin_ID, by = "VarID")
####################################################
######## create data set with all info #############
bp_windows = list("5", "10", "15", "20")
####################################################
############ make an output dataframe ##############
output = bin_indels 
####################################################
############### calculate the rAF ##################
for (bp in bp_windows) {
  ## first step is to find the mean AN of the indels that share the same CHR number and groupID number for a certain bp window 
  ## whenever bp is in that window, create a new data frame called mean_AN that keeps track of the mean_AN
  if (bp == "5") {
    ## aggregate by CHR and GROUPID num for bp 5
    mean_AN = aggregate(bin_indels$AN, by = list(bin_indels$GROUPID_bp5, bin_indels$CHR), mean)
  }
  if (bp == "10") {
    ## aggregate by CHR and GROUPID num for bp 10
    mean_AN = aggregate(bin_indels$AN, by = list(bin_indels$GROUPID_bp10, bin_indels$CHR), mean)
  }
  if (bp == "15") {
    ## aggregate by CHR and GROUPID num for bp 15
    mean_AN = aggregate(bin_indels$AN, by = list(bin_indels$GROUPID_bp15, bin_indels$CHR), mean)
  }
  if (bp == "20") {
    ## aggregate by CHR and GROUPID num for bp 20
    mean_AN = aggregate(bin_indels$AN, by = list(bin_indels$GROUPID_bp20, bin_indels$CHR), mean)
  }
  
  ## set variable that will change as the bp changes 
  groupID_bpx = paste0("GROUPID_bp",bp)
  meanAN_bpx = paste0("meanAN_", bp)
  
  ## rename the meanAN dataframe
  colnames(mean_AN) = c(groupID_bpx, "CHR", meanAN_bpx)
  mean_AN = mean_AN[, c("CHR", groupID_bpx, meanAN_bpx)]
  
  
  ## merge with the indels dataframe by CHR and GROUPID for specific bp to add the meanAN column into the big dataframe
  output = merge(output, mean_AN, by=c("CHR", groupID_bpx))
  ## for gnomAD data: rAF is 2(allelecount/meanAN) so the threshold is 1 * 10^-4 
  ## finally calculate the rAF by dividing the AC of the bp window by the mean AN for the bp window 
  if (database == "gnomAD") {
    if (bp == "5") {
      output$rAF5 = 2*(output$AlleleCount_bp5 /output$meanAN_5)
    }
    if (bp == "10") {
      output$rAF10 = 2*(output$AlleleCount_bp10 /output$meanAN_10)
    }
    if (bp == "15") {
      output$rAF15 = 2*(output$AlleleCount_bp15 /output$meanAN_15)
    }
    if (bp == "20") {
      output$rAF20 = 2*(output$AlleleCount_bp20 /output$meanAN_20)
    }
  }
  
  ## for IGM data: rAF is (allelecount/meanAN) so the threshold is 1 * 10^-4 
  if (database == "IGM") {
    if (bp == "5") {
      output$rAF5 = (output$AlleleCount_bp5 /output$meanAN_5)
    }
    if (bp == "10") {
      output$rAF10 = (output$AlleleCount_bp10 /output$meanAN_10)
    }
    if (bp == "15") {
      output$rAF15 = (output$AlleleCount_bp15 /output$meanAN_15)
    }
    if (bp == "20") {
      output$rAF20 = (output$AlleleCount_bp20 /output$meanAN_20)
    }
  }
}

## sanity check:
## compare output with indels_test 
indels_test$VarID = paste0(indels_test$CHR, "-", indels_test$POS, "-", indels_test$REF, "-", indels_test$ALT)
indels_test = indels_test[order(indels_test$VarID),]
output = output[order(output$VarID),]

indels_test$raf5_diff = abs(indels_test$rAF5 - output$rAF5)
indels_test$raf10_diff = abs(indels_test$rAF10 - output$rAF10)
indels_test$raf15_diff = abs(indels_test$rAF15 - output$rAF15)
indels_test$raf20_diff = abs(indels_test$rAF20 - output$rAF20)

## The biggest difference is 10^-15, which is likely due to rounding errors. 

## format CSV output 
output = output %>% mutate (
  CHR = CHR, 
  POS = POS, 
  REF = REF, 
  ALT = ALT, 
  AC = AC, 
  GROUPID_bp5 = GROUPID_bp5, 
  GROUPID_bp10 = GROUPID_bp10, 
  GROUPID_bp15 = GROUPID_bp15,
  GROUPID_bp20 = GROUPID_bp20, 
  AlleleCount_bp5 = AlleleCount_bp5,
  AlleleCount_bp10 = AlleleCount_bp10,
  AlleleCount_bp15 = AlleleCount_bp15, 
  AlleleCount_bp20 = AlleleCount_bp20, 
  AN = AN, 
  sAF = sAF,
  AN5 = meanAN_5,
  AN10 = meanAN_10,
  AN15 = meanAN_15,
  AN20 = meanAN_20,
  VarID = VarID,
  rAF5 = rAF5,
  rAF10 = rAF10,
  rAF15 = rAF15,
  rAF20 = rAF20
) %>% select(CHR, POS, REF, ALT, AC, GROUPID_bp5, GROUPID_bp10, GROUPID_bp15, GROUPID_bp20, AlleleCount_bp5, AlleleCount_bp10, AlleleCount_bp15, AlleleCount_bp20, AN, sAF, AN5, AN10, AN15, AN20, VarID, rAF5, rAF10, rAF15, rAF20)
## write CSV output 

fwrite(output, paste0( Sys.Date(), "_", database, "_sAF_rAF.csv"), row.names = FALSE, quote = FALSE)












