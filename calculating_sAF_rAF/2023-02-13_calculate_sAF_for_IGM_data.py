## Sandy Yang
## September 6 2022
## This is the code I ran to calculate the sAF for the IGM data

import pandas as pd
import optparse
import csv 

## commandline inputs
parser = optparse.OptionParser()
parser.add_option('-a', '--AN_file', action='store', dest="AN_file")
options, args = parser.parse_args()

AN_file = pd.read_csv(options.AN_file, sep=',')
AN_file.columns = ['CHR','POS','REF','ALT','AC','GROUPID_bp5','GROUPID_bp10','GROUPID_bp15','GROUPID_bp20','AlleleCount_bp5','AlleleCount_bp10','AlleleCount_bp15','AlleleCount_bp20','VarID','AN']

fields = ['CHR', 'POS', 'REF', 'ALT', 'AC', 'AN', 'sAF']

## Create an output .csv file that separates the VarID into its components while also including a sAC column 
with open ('IGM_indels_with_sAF.csv', 'w') as output:
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
		## calculate the sAF
		## Note: for IGM, we are finding the sAF by multiply the AN by 2. This is equivalent to dividing the whole sAF value by 2. 
		sAF = (1/2)*(int(AN_file['AC'][i])/int(AN_file['AN'][i]))

		variant = []

		variant.append(AN_file['CHR'][i])
		variant.append(AN_file['POS'][i])
		variant.append(AN_file['REF'][i])
		variant.append(AN_file['ALT'][i])
		variant.append(AN_file['AC'][i])
		variant.append(AN_file['AN'][i])
		variant.append(sAF)

		row.append(variant)
	writer.writerows(row)

## command: python3 2023-02-13_calculate_sAF_for_IGM_data.py -a 20220906_IGM__VarID_AC_AN_no_duplic_sorted_header.csv
