# Introduction to regional allele frequency (rAF)
To increase the reliability of allele frequency values of insertion-deletion polymorphisms (indels) in genetic diagnostics and gene discovery efforts, we proposed the use of regional allele frequency (rAF) in lieu of population allele frequency (AF or sAF). Regional allele frequency takes into account the number of indels within a base-pair window (10bp, 20bp, 30bp, 40bp) in order to offset errors in AF estimations from differences in mapping and variant-calling methods. 

# Calculation of standard and regional allele frequencies
Population allele frequency (also referred to as "standard allele frequency" or sAF) is defined by the allele count (AC) of a particular variant, divided by the number of covered alleles (AN). 

The regional allele frequency (rAF) is defined by the regional allele count of a particular variant, divided by the regional AN. The regional allele count is the number of variants within a particular base-pair window (10bp, 20bp, 30bp, 40bp) and can be calculated by finding the sum of the AC values of the indels within the window. The regional AN is the mean of the AN values of the variants within a particular base-pair window (10bp, 20bp, 30bp, 40bp).

# Discrepancy between gnomAD AN values and IGM AN values
The AN values of the IGM data only accounts for the number of individuals. However, the AN values of the gnomAD data represent the number of covered alternate alleles. Therefore, it is crucial to remember to multiply the IGM AN number by 2 during calculations of allele frequency. 

# Rareness
Rare is defined as allele frequency that is less than or equal to (<=) 10^-4 for both gnomAD and IGM data. 

Common is defined as allele frequency that is greater than (>) 10^-4 for both gnomAD and IGM data. 

# The Genome Aggregation Database (gnomAD)
Source: https://gnomad.broadinstitute.org
Version: v2.1.1

# Tutorial 
<a href="https://github.com/ColumbiaCPMG/RegionalAlleleFrequency/blob/main/tutorial.md" target="_blank">Begin Here</a>
