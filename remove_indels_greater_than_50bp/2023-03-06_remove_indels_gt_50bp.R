## Sandy Yang
## March 6, 2023 
## This script removes indels that are > 50 bp and only keeps indels <= 50bp 

## Read in files 
gnomad_indels = read.csv("2023-02-17_gnomAD_sAF_rAF.csv", sep=',', header = TRUE)
IGM_indels = read.csv("finalized_project_data/2023-02-17_IGM_sAF_rAF.csv", sep=',', header = TRUE)

## remove the indels that have > 50bp and only keep indels <= 50bp in length 
gnomad_indels_lte_50 = gnomad_indels %>% mutate(indel_length = nchar(gnomad_indels$ALT) - nchar(gnomad_indels$REF)) %>% filter(indel_length <= 50)
IGM_indels_lte_50 = IGM_indels %>% mutate(indel_length = nchar(IGM_indels$ALT) - nchar(IGM_indels$REF)) %>% filter(indel_length <= 50)

## Write files 
fwrite(gnomad_indels_lte_50, paste0( Sys.Date(), "_gnomAD_sAF_rAF_lte_50bp.csv"), row.names = FALSE, quote = FALSE)

fwrite(IGM_indels_lte_50, paste0( Sys.Date(), "_IGM_sAF_rAF_lte_50bp.csv"), row.names = FALSE, quote = FALSE)
