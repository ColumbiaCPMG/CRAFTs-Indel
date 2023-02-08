import numpy as np
import pandas as pd 
import seaborn as sns
import pickle
import matplotlib
import matplotlib.pyplot as plt
import statsmodels
df = pd.read_csv("2022-12-30_figure5_gnomad_input.csv")
filename = 'gnomAD' #or IGM
#bp_window = '5'  # options: 5,10,15,20 We already know ours is bp window 5 because that is our input 
#log_x = np.asarray(log_df['rAF'+str(bp_window)])
log_x = df['rAF5']
log_y = df['log_key']
#log_y = np.asarray(log_df['log_key'])
sns.regplot(x = log_x, y = log_y, logistic = True, color = "gray")
plt.xlim([0,1])
plt.xlabel('Regional Allele Frequency of 10 bps Window', fontsize = 20)
plt.ylabel('Log Odds (Pathogenic/Likely Pathogenic)', fontsize =20)
plt.title(str(filename)+': The logistic regression curve of benign/ pathogenic indels (10bps window) \n Pathogenic/Likely Pathogenic = 1, Benign/Likely Benign = 0', fontsize=20 )
plt.tick_params(axis='both', which='major', labelsize=20)
plt.show()
