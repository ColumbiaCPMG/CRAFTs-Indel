import pickle
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#bp_window = '5'  # options: 5,10,15,20
filename = 'gnomAD' #or IGM
## Read in csv file with Clinvar labels
df = pd.read_csv("2022-12-30_figure5_gnomad_input.csv")
clinvar_class = ['B/LB','P/LP']

## first, to make our graph, we only want to look at rAF5 < 0.05 (a zoomed up view of our graph)
df_0_05 = df.loc[df['rAF5'] <= 0.05]

fig, ax = plt.subplots()
rects1 = ax.hist([df_0_05.loc[df_0_05['clinvar_class'] == 'B/LB']['rAF5'],df_0_05.loc[df_0_05['clinvar_class'] == 'P/LP']['rAF5']],
                 bins =20, log = True,label = ['Benign/Likely Benign','Pathogenic/Likely Pathogenic'], color = ["gray","black"])
ax.set_ylabel('Number of Indels on a Log Scale', fontsize = 20)
ax.set_xlabel('Regional Allele Frequency (Window = 10bps)', fontsize = 20)
ax.set_title(str(filename)+': The rAF distribution of pathogenic/benign indels (window = 10bps)', fontsize = 20)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize = 20)
plt.tick_params(axis='both', which='major', labelsize=20)
fig.tight_layout()
plt.show()

