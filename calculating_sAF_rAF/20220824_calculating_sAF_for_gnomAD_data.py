## Author: Sandy Yang
## Date: August 24 2022
## This code will be part of the effort to calculate the rAF and the sAF given the rAC and the sAC 

## To calculate the sAF, this would be the sAC/AN 
## To calcuate the rAF, this would be the rAF/(mean AN for that bin)

## Here begins the code to calculate the sAF 
## Read in the the 2022_08_23_gnomAD2.1.1_VarIDs_AC_AN_Indels.csv

import pandas as pd
import optparse
import csv 

## commandline inputs
parser = optparse.OptionParser()
parser.add_option('-a', '--AN_file', action='store', dest="AN_file")
options, args = parser.parse_args()

AN_file = pd.read_csv(options.AN_file, sep=',')
AN_file.columns = ['VarID', 'AC', 'AN', 'TYPE']

AN_file = AN_file[AN_file['AN'] != 0]
AN_file = AN_file[AN_file['AC'] != 0]
AN_file = AN_file.reset_index(drop= True)

fields = ['CHR', 'POS', 'REF', 'ALT', 'AC', 'AN', 'sAF']

## Create an output .csv file that separates the VarID into its components while also including a sAC column

with open ('gnomad_indels_with_sAF.csv', 'w') as output:
	## creating csv writer object
	writer = csv.writer(output)
	## write the column fields
	writer.writerow(fields)
	## initialize empty list 
	row = []
	## for every row in the input csv file with the AN information
	for i in range(0, len(AN_file['VarID'])):
	## for testing: 
	## for i in range(0,100):
		## first, separate the VarID field into its components
		VarID_split = AN_file['VarID'][i].split('\t')
		##print (VarID_split)
		CHR = VarID_split[0]
		POS = VarID_split[1]
		REF = VarID_split[2]
		ALT = VarID_split[3]

		## calculate the sAF
		## Note: for gnomAD, we are finding the sAF by dividing the AN by 2 (equivalent to multiplying the AC/AN by 2) 
		sAF = 2*(AN_file['AC'][i]/AN_file['AN'][i])

		variant = []

		variant.append(CHR)
		variant.append(POS)
		variant.append(REF)
		variant.append(ALT)
		variant.append(AN_file['AC'][i])
		variant.append(AN_file['AN'][i])
		variant.append(sAF)

		row.append(variant)
	writer.writerows(row)










