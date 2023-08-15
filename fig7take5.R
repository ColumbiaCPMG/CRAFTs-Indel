#Figure 7 take 4
install.packages("circlize")
library(circlize)

install.packages("BiocManager")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
library(GenomicRanges)

install.packages("data.table")
library(data.table)


setwd("/Users/User_1/Desktop/rAF_project_2")

#####load in data#####

##rAF-hi 10bp 
rAF_hi_bp_10_df_1 = fread("gnomad.exomes.r2.1.1.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_10_df_1$Chr <- paste("chr", rAF_hi_bp_10_df_1$Chr, sep = "")
colnames(rAF_hi_bp_10_df_1) <- c("chr", "start", "end")

#rAF_hi_bp_10_df_1$start_end <- paste(rAF_hi_bp_10_df_1$start, "-", rAF_hi_bp_10_df_1$end, sep = "")


#testing if this will work when start and end pos are in chr format
# Convert a column from integer to character
# rAF_hi_bp_10_df_1$start <- as.character(rAF_hi_bp_10_df_1$start)
# rAF_hi_bp_10_df_1$end <- as.character(rAF_hi_bp_10_df_1$end)




rAF_hi_bp_10_df_2 = fread("2023-03-23_IGM_n39367_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_10_df_2$Chr <- paste("chr", rAF_hi_bp_10_df_2$Chr, sep = "")
colnames(rAF_hi_bp_10_df_2) <- c("chr", "start", "end")

#UK.BB liftover rAF-hi 10bp 
rAF_hi_bp_10_df_3 = fread("UK.BB.exomes.430k.sites_indelsonly_rAF_bp10_rAF_hiIndels.lt50bp.region.bed")
rAF_hi_bp_10_df_3 = rAF_hi_bp_10_df_3[, c("V1", "V2", "V3")]
colnames(rAF_hi_bp_10_df_3) <- c("chr", "start", "end")

#LCR-file 
lcr_file_37=fread('numInsteadOfGRCh37_AllTandemRepeatsandHomopolymers_slop5.bed')



####################################################
####################################################
####################################################
####################################################



bed_data <- rAF_hi_bp_10_df_1
gr <- GRanges(seqnames = bed_data$chr, ranges = IRanges(start = bed_data$start, end = bed_data$end))
circos.initializeWithIdeogram()

circos.genomicTrack(
  bed_data,
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(gr, col = "blue")
  }
)




bed_data <- read.table("path/to/your/file.bed", header = FALSE)

# Set up Circos parameters
circos.par(cell.padding = c(0, 0, 0, 0))

# Create a Circos plot
circos.genomicPlot(
  bed_data,
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(region, col = "blue")
  }
)











circos.par(cell.padding = c(0, 0, 0, 0))
circos.initialize(sectors = bed_data$chr)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {})
circos.clear()


# Example heatmap track
circos.genomicHeatmap(gr, col = c("blue", "white", "red"), width = 0.05)

# Example link track
circos.link(gr, gr, col = "blue")

# Example histogram track
circos.histogram(gr)


# Save the plot to a file (e.g., PNG)
png("circos_plot.png", width = 800, height = 800)
circos.clear()
dev.off()


circos.clear()

