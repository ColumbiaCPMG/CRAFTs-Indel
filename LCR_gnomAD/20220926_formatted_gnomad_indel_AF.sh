## January 6, 2023 
## Sandy Yang
## This code formats the input file 
## CHR,POS,VarID,sAF,rAF5,rAF10,rAF15,rAF20

awk -F "," '{if (NR==1) {
	print $1","$2","$20","$15","$21","$22","$23","$24
}
if (NR!=1) {
	print $1","$2","$20","$15","$21","$22","$23","$24}
}' 2023-03-06_gnomAD_sAF_rAF_lte_50bp.csv > 2023-03-06_gnomad_indels_AF.csv