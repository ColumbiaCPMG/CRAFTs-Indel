library(tidyverse)
library(data.table)
library(readxl)
library(dplyr)
library(stringr)

#install.packages("circlize")
library(circlize)

# Load your BED data 

setwd ("")

############################################################
## Convert UKBB 10bp rAF-hi regions from GRCh38 to GRCh37 
############################################################
UKBB_10bp_rAF_hi_38 = read.table("UKBB/UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.region.bed", header = FALSE)
UKBB_10bp_rAF_hi_38$V1 = paste0("chr",UKBB_10bp_rAF_hi_38$V1)
fwrite(UKBB_10bp_rAF_hi_38, "UKBB_10bp_rAF_hi_GRCh38.bed", sep = "\t", col.names = F)

## Use this website: https://genome.ucsc.edu/cgi-bin/hgLiftOver
## min number of bases to revamp is 0.1
## liftover from grch38 to grch37
## save it in a new file and that will be your UKBB input 

############################################################
## Read in appropriate bedfiles for LCR, gnomAD, IGM, UKBB 
############################################################

LCR_GRCh37 = read.table("LCR/GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed", header = FALSE)
gnomAD_10bp_rAF_hi = read.table("gnomAD/gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.region.bed", header = T)
IGM_10bp_rAF_hi = read.table("IGM/2023-03-23_IGM_n39367_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.region.bed", header = T)
UKBB_10bp_rAF_hi = read.table("UKBB/UKBB_rAF_hi_regions_10bp_CHGR37.bed" , header = FALSE)

############################################################
## Format bedfiles
############################################################

## Format LCR file 
LCR_GRCh37$chr = paste0("chr",LCR_GRCh37$V1)
LCR_GRCh37$start = LCR_GRCh37$V2
LCR_GRCh37$end = LCR_GRCh37$V3
LCR_GRCh37$numeric_column = 1
LCR_GRCh37 = LCR_GRCh37 %>% select("chr", "start", "end", "numeric_column")

## Format gnomAD bedfile 
gnomAD_10bp_rAF_hi$Chr = paste0("chr", gnomAD_10bp_rAF_hi$Chr)
gnomAD_10bp_rAF_hi$numeric_column = 1
colnames(gnomAD_10bp_rAF_hi) = c("chr", "start", "end", "numeric_column")

## Format IGM bedfile 
IGM_10bp_rAF_hi$Chr = paste0("chr", IGM_10bp_rAF_hi$Chr)
IGM_10bp_rAF_hi$numeric_column = 1
colnames(IGM_10bp_rAF_hi) = c("chr", "start", "end", "numeric_column")

## Format UKBB bedfile 
UKBB_10bp_rAF_hi$numeric_column = 1
UKBB_10bp_rAF_hi = UKBB_10bp_rAF_hi %>% select("V1", "V2", "V3", "numeric_column")
colnames(UKBB_10bp_rAF_hi) = c("chr", "start", "end", "numeric_columns")


############################################################
## create circle plot 
############################################################
 
circos.clear()

## Note: by default, circos.initializeWithIdeogram() uses Hg19/GRCh37 
## initialize the ideogram in Hg19 
circos.initializeWithIdeogram(species = "hg19")
circos.par("track.height" = 0.1)

## label ideogram
#circos.labels("chr18", 1, "A", side = "inside")

## create a track with the LCR for Hg19 
circos.genomicTrack(LCR_GRCh37, ylim= c(0,1), numeric.column = 4,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h")
                    })
## label LCR ring
#circos.labels("chr19", 1, "B", side = "inside")
print("LCR finished")


## create a track with the gnomAD rAF-hi regions by 10bp window 
circos.genomicTrack(gnomAD_10bp_rAF_hi, ylim= c(0,1), numeric.column = 4,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h")
                    })
## label gnomAD ring
#circos.labels("chr20", 1, "C", side = "inside")
print("gnomAD finished")

## create a track with the IGM rAF-hi regions by 10bp window 
circos.genomicTrack(IGM_10bp_rAF_hi, ylim= c(0,1), numeric.column = 4,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h")
                    })
## label IGM ring 
#circos.labels("chr21", 1, "D", side = "inside")
print("IGM finished")

## create a track with the UKBB rAF-hi regions by 10bp window 
circos.genomicTrack(UKBB_10bp_rAF_hi, ylim= c(0,1), numeric.column = 4,
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h")
                    })
## label UKBB ring 
#circos.labels("chr22", 1, "E", side = "inside")
print("UKBB finished")

pdf("fig7.pdf")
circos.clear()
dev.off()

############################################################
## Make a key 
############################################################

#circos.clear()

#sectors = unique(gnomAD_10bp_rAF_hi$chr)
#sectors = c("chr1")

#circos.initializeWithIdeogram(species = "hg19")
#circos.par("track.height" = 0.1)

#circos.track(sectors, ylim = c(0,1))
#circos.text(1, 0.5, c("LCR"), niceFacing =T)

#circos.track(sectors, ylim = c(0,1))
#circos.text(1, 0.5, c("gnomAD"), niceFacing =T)

#circos.track(sectors, ylim = c(0,1))
#circos.text(1, 0.5, c("IGM"), niceFacing =T)

#circos.track(sectors, ylim = c(0,1))
#circos.text(1, 0.5, c("UKBB"), niceFacing =T)

#circos.clear()
#dev.off()

