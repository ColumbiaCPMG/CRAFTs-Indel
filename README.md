# Paper
<a href="https://link.springer.com/epdf/10.1007/s10142-024-01358-3?sharing_token=fOQMmH39idz4_03JoizO9_e4RwlQNchNByi7wbcMAY7f5dw06VKtYPRgIiI2hEdIVysBRgLtOXQj7PPVNBpV9XhVjDqC4X0tew3FuhV8b4H1JpL76P-OQKBVHujeQyrAO6wYZnVMA7sGMKyoIsZhhPiWIELt9URtLTWzkd8llcU%3D" target="_blank">Assisting theanalysis ofinsertions anddeletions using regional allele frequencies</a>

# Contributors
Sarath Babu Krishna Murthy, Sandy Yang, Shiraz Bheda, Nikita Tomar, Haiyue Li, Amir Yaghoobi, Atlas Khan, Krzysztof Kiryluk, Joshua E. Motelow, Nick Ren, Ali G. Gharavi, Hila Milo Rasouly.
 
# Introduction to regional allele frequency (rAF)
To increase the reliability of allele frequency values of insertion-deletion polymorphisms (indels) in genetic diagnostics and gene discovery efforts, we proposed the use of regional allele frequency (rAF) in lieu of population allele frequency (AF or sAF). Regional allele frequency takes into account the number of indels within a base-pair window (10bp, 20bp, 30bp, 40bp) in order to offset errors in AF estimations from differences in mapping and variant-calling methods. 

# Calculation of standard and regional allele frequencies
Population allele frequency (also referred to as "standard allele frequency" or sAF) is defined by the allele count (AC) of a particular variant, divided by the number of covered alleles (AN). 

The regional allele frequency (rAF) is defined by the regional allele count of a particular variant, divided by the regional AN. The regional allele count is the number of variants within a particular base-pair window (10bp, 20bp, 30bp, 40bp) and can be calculated by finding the sum of the AC values of the indels within the window. The regional AN is the mean of the AN values of the variants within a particular base-pair window (10bp, 20bp, 30bp, 40bp).

# Note on IGM AN values
ATAV which holds the IGM dataset provides the number of sampels in which a particular variant is covered. So in order to use the AN value for sAF/rAF calculation the number samples covered for a variant is multiplied by 2 accounting for both alleles. 

# Rareness
In this study, Rare is defined as allele frequency that is less than or equal to (<=) 10^-4 and Common is defined as allele frequency that is greater than (>) 10^-4. This is calculated based on the size of the cohort. This value of rareness should be adjusted according to study. 


# The Genome Aggregation Database (gnomAD)
Source: https://gnomad.broadinstitute.org
Version: v2.1.1

# UK Biobank
Source: https://www.ukbiobank.ac.uk/

# Index 
<a href="https://github.com/ColumbiaCPMG/RegionalAlleleFrequency/blob/main/tutorial.md" target="_blank">Begin Here</a>
