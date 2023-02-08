## January 6, 2023 
## Sandy Yang
## This code formats the input file 
## CHR,POS,VarID,sAF,rAF5,rAF10,rAF15,rAF20

awk -F "," '{if (NR==1) {
	print $1","$2","$20","$15","$21","$22","$23","$24
}
if (NR!=1) {
	print $1","$2","$20","$15","$21","$22","$23","$24}
}' 2022-12-30_gnomAD_sAF_rAF.csv  > 2023-01-06_gnomad_indels_AF.csv