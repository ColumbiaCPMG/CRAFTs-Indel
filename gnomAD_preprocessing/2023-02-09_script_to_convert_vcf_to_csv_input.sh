## code to get the csv file from the vc file downloaded from gnomad website 

bcftools view --include 'TYPE="INDEL"' gnomad.exomes.r2.1.1.sites.vcf.bgz -O b >gnomad.exomes.r2.1.1.sites_indelsonly.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AC\n' gnomad.exomes.r2.1.1.sites_indelsonly.vcf.gz >gnomad.exomes.r2.1.1.sites_indelsonly_AC.txt

## code to get the AN, AC columns 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT,%AC,%AN,%TYPE\n' gnomad.exomes.r2.1.1.sites_indelsonly.vcf.gz >gnomAD2.1.1_VarIDs_AC_AN_Indels.csv

 ( echo "VarID,AC,AN,TYPE" ; cat gnomAD2.1.1_VarIDs_AC_AN_Indels.csv ) > gnomAD2.1.1_VarIDs_AC_AN_Indels_header.csv