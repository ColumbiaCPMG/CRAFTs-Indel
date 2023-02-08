## Author: Haiyue "Arianna" Li
## This code is to find the pickle file for Fig 1
## This code also calculates bin length, which is needed for Fig 1 

# variables
filename = 'gnomAD' #the name of the data working with
rare = 1 * 10 ** -4  #rare threshold of your choice
bps = [5, 10, 15, 20] #choose the rAF bp, note this is half of the sliding window size
verbose = True #if you wish to print out intermediate steps 
long = 140 #min amount of bps spanned by bin for a bin to be classifed as long bin

# load packages
import pandas as pd
import numpy as np
from collections import Counter
pd.set_option('display.max_columns', None)

#load data (everything) 
data = pd.read_csv(r"2022-12-30_gnomAD_sAF_rAF.csv",
            sep = ",",na_filter = False,low_memory = False)

print("Finished loading in data.")

## This part can take a long time for the gnomAD data (> 2 hours)
# Add indel type (ins or del) and insertion/deletion lens to the end of file data
data['type'] = None
data['ins_len'] = None
for index, row in data.iterrows():
    data.at[index, 'ins_len'] = len(data['ALT'][index]) - len(data['REF'][index])
    if len(data['ALT'][index]) > len(data['REF'][index]):
        data.at[index, "type"] = "insertion"
    else:
        data.at[index, "type"] = "deletion"

print("Finished determining indel type (insertion/deletion).")


# To create data_50, indels whose lengths >=50bps are removed
# data_50 have indels that have less than 50bp in length 
data_50 = data[data['ins_len'] <= 50]
data_50 = data_50[data_50['ins_len'] >= -50]


# separate data_50 into chromosomes to get data_c50
## determine a list of all the chromosomes and order the X, Y and MT chromosomes 
chr_unique = np.unique(data['CHR'])
c_list = chr_unique.tolist()
c_list = [item.replace('X', '98') for item in c_list]
c_list = [item.replace('Y', '99') for item in c_list]
c_list = [item.replace('MT', '97') for item in c_list]
c_list.sort(key=int)
c_list = [item.replace('98', 'X') for item in c_list]
c_list = [item.replace('99', 'Y') for item in c_list]
c_list = [item.replace('97', 'MT') for item in c_list]

## data_c50 contains all the indels that are less than 50bp in length and are organized in a dictionary by their chr number 
data_c50 = {}
for c in c_list:
    data_c50[c] = pd.DataFrame()
    data_c50[c] = data_50.loc[data_50['CHR'] == c]

data_c50 = {k: data_c50[k] for k in c_list}

for c in c_list:
    data_c50[c] = data_c50[c].reset_index()

## We are adding a col to determine whether they are rare by sAF/rAF
## These columns are empty columns currently
for c in c_list:
    for bp in bps:
        data_c50[c]['standard rare'] = None
        data_c50[c]['regional rare' + str(bp)] = None

# Determine if an indel is rare or common by standard AF and populate the empty columns 
for c in c_list:
    if verbose:
        print(c)
    for index, row in data_c50[c].iterrows():
        if data_c50[c].at[index, "sAF"] <= rare:
            data_c50[c].at[index, "standard rare"] = 'rare'
        else:
            data_c50[c].at[index, "standard rare"] = 'com'

# Determine if an indel is rare or common by regional AF and populate the empty columns 
for c in c_list:
    for bp in bps:
        if verbose:
            print("chr:" + str(c) + " bp:" + str(bp))
        for index, row in data_c50[c].iterrows():
            if data_c50[c].at[index, "rAF" + str(bp)] <= rare:
                data_c50[c].at[index, "regional rare" + str(bp)] = 'rare'
            else:
                data_c50[c].at[index, "regional rare" + str(bp)] = 'com'

## data_z50 are indels that are rare sAF/common rAF only and also less than 50bp in length 
data_z50 = {}
for bp in bps:
    data_z50[bp] = {}
    for c in c_list:
        data_z50[bp][c] = data_c50[c][
            (data_c50[c]["regional rare" + str(bp)] == 'com') & (data_c50[c]["standard rare"] == 'rare')]

z50_binid = {}
for bp in bps:
    z50_binid[bp] = {}
    for c in c_list:
        z50_binid[bp][c] = np.unique(data_z50[bp][c]['GROUPID_bp' + str(bp)]).tolist()


# data_a50 is to find the indels that are in a bin with at least one rare sAF/common rAF indel 
# all indels (not just rare sAF; common rAF)
# data_a50:
## the bin has to contain at least one indel that is rare sAF/common rAF, all the indels in this bin are kept regardless of whether they are also rare sAF/common rAF as long as one indel in that bin is 
## we are keeping the indels that are in bins with at least one indel that is rare sAF/common rAF because this helps to maintain the integrity of the bin
## if we removed these indels, then indels that are in a bin will no longer be close to each other because the position where you begin to count from will change if you disregard these indels 

data_a50 = {}
for bp in bps:
    data_a50[bp] = {}
    for c in c_list:
        data_a50[bp][c] = data_c50[c][data_c50[c]["GROUPID_bp" + str(bp)].isin(z50_binid[bp][c])]

for bp in bps:
    for c in c_list:
        data_a50[bp][c]['interest'] = None

for bp in bps:
    for c in c_list:
        if verbose:
            print("chr:" + str(c) + " bp:" + str(bp))
        for index, row in data_a50[bp][c].iterrows():
            if data_a50[bp][c].at[index, "standard rare"] == 'rare' and data_a50[bp][c].at[
                index, "regional rare" + str(bp)]:
                data_a50[bp][c].at[index, "interest"] = 'yes'
            else:
                data_a50[bp][c].at[index, "interest"] = 'no'


"""# 3 find data_a50e: same as data_a, but remove bins containing only 1 indel"""

data_a50e = {}  # e for edited, bins containing only 1 indel is removed
for bp in bps:
    data_a50e[bp] = {}
    for c in c_list:
        data_a50e[bp][c] = pd.DataFrame()

for bp in bps:
    for c in c_list:
        data_a50[bp][c] = data_a50[bp][c].sort_values(by=['POS'])

## sort by group ID 
for c in c_list:
    data_a50e[5][c] = data_a50[5][c].sort_values(by=['GROUPID_bp5'])
    data_a50e[5][c] = data_a50[5][c][data_a50[5][c].GROUPID_bp5.duplicated(keep=False)]
    data_a50e[10][c] = data_a50[10][c].sort_values(by=['GROUPID_bp10'])
    data_a50e[10][c] = data_a50[10][c][data_a50[10][c].GROUPID_bp10.duplicated(keep=False)]
    data_a50e[15][c] = data_a50[15][c].sort_values(by=['GROUPID_bp15'])
    data_a50e[15][c] = data_a50[15][c][data_a50[15][c].GROUPID_bp15.duplicated(keep=False)]
    data_a50e[20][c] = data_a50[20][c].sort_values(by=['GROUPID_bp20'])
    data_a50e[20][c] = data_a50[20][c][data_a50[20][c].GROUPID_bp20.duplicated(keep=False)]

"""# 5 long bins: based on bp count per bin"""

''' for data_a, some bins are left with just 1 indel, so calculate only data_a50e instead'''

# create empty bins to calculate bin length (bp count of bin) 
data_bin = {}
for bp in bps:
    data_bin[bp] = {}
    for c in c_list:
        data_bin[bp][c] = pd.DataFrame()

# output: data_bin[c][p] - an overall dataframe that contains all info for start and end indels of each bin
## within each chromosome, 
## sort values first by position 
for c in c_list:
    for bp in bps:
        ### should this be data_bin[bp][c], not data_a50e[bp][c]????
        data_a50e[bp][c] = data_a50e[bp][c].sort_values(by='POS')
## group bins by groupID and set to data_bin[bp][c]
for c in c_list:
    for bp in bps:
        data_bin[bp][c] = data_a50e[bp][c].groupby(by='GROUPID_bp' + str(bp), as_index=False).nth([0, -1])

## within the bins, sort the groupID again (same bin has same groupID) 
for bp in bps:
    for c in c_list:
        data_bin[bp][c] = data_bin[bp][c].sort_values(by="GROUPID_bp" + str(bp))
## within the bins, sort by POS 
for c in c_list:
    for bp in bps:
        data_bin[bp][c] = data_bin[bp][c].sort_values(by="POS")
        data_bin[bp][c] = data_bin[bp][c].reset_index()

## create an empty col called condition to identify which indel is starting/ending indel 
for c in c_list:
    for bp in bps:
        data_bin[bp][c]['condition'] = None

## Look for start and end of the indels within the bins 
# for each chromosome 
for c in c_list:
    # for each bp 
    for bp in bps:
        # each row in data_bin[bp][c] would be all the indels grouped by (aggregated by) groupID and sorted by groupID and pos 
        for index, row in data_bin[bp][c].iterrows():
            # if the index is even, then the indel is the starting indel?
            if index % 2 == 0:
                data_bin[bp][c].at[index, "condition"] = "start"
            else:
                data_bin[bp][c].at[index, "condition"] = "end"

print("Beginning bin length calculation.")

# bin length calculation
# create empty column for bin length
for c in c_list:
    for bp in bps:
        data_bin[bp][c]['bin_len'] = None

# calculate and populate bin length column 
for c in c_list:
    for bp in bps:
        if verbose:
            print("chr:" + str(c) + " bp:" + str(bp))
        for index, row in data_bin[bp][c].iterrows():
            a = data_bin[bp][c][['POS']].diff()
            data_bin[bp][c]["bin_len"] = a

for c in c_list:
    for bp in bps:
        print("chr:" + str(c) + " bp:" + str(bp))
        for index, row in data_bin[bp][c].iterrows():
            if data_bin[bp][c].at[index, "condition"] == "end":
                if data_bin[bp][c].at[index, "type"] == "insertion":
                    data_bin[bp][c].at[index, "bin_len"] = data_bin[bp][c].at[index, "bin_len"] + data_bin[bp][c].at[
                        index, "ins_len"]

for c in c_list:
    for bp in bps:
        for index, row in data_bin[bp][c].iterrows():
            if data_bin[bp][c].at[index, "condition"] == "start":
                data_bin[bp][c].at[index, "bin_len"] = np.NaN

for c in c_list:
    for bp in bps:
        data_bin[bp][c] = data_bin[bp][c].fillna(method="backfill", axis=0)

print("Finished bin length calculation.")


print("Beginning to export input file for Figure 1 script.")


# count how many bins there are for bins with different bin length
import pickle
a3 = {}
for bp in bps:
    a3[bp] = []
    for c in c_list:
        a1 = data_bin[bp][c][data_bin[bp][c]['condition'] == 'start']
        a2 = a1['bin_len'].tolist()
        a3[bp] = a3[bp] + a2
handle = open('2022-12-07_gnomad_fig1_input.pkl', 'wb')
pickle.dump(a3, handle)

print("Finished exporting input file for Figure 1 script.")
