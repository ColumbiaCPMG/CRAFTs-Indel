## January 6, 2023 
## In an effort to try to make the LCR analysis go quicker, I'm going to try to use bedtools intersect 


## For only rare by standard and common by regional 
##bp 5 
bedtools intersect -a 2023-03-06_rare_sAF_common_rAF_5_bedfile.bed -b GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed | sort -k1,1n -k2,2n > 2023-03-06_LCR_rare_sAF_common_rAF_5.bed

##bp 10
bedtools intersect -a 2023-03-06_rare_sAF_common_rAF_10_bedfile.bed -b GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed | sort -k1,1n -k2,2n > 2023-03-06_LCR_rare_sAF_common_rAF_10.bed

##bp 15 
bedtools intersect -a 2023-03-06_rare_sAF_common_rAF_15_bedfile.bed -b GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed | sort -k1,1n -k2,2n > 2023-03-06_LCR_rare_sAF_common_rAF_15.bed

##bp 20
bedtools intersect -a 2023-03-06_rare_sAF_common_rAF_20_bedfile.bed -b GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed | sort -k1,1n -k2,2n > 2023-03-06_LCR_rare_sAF_common_rAF_20.bed