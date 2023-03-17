import pickle
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

filename = 'gnomAD' 
## Read in csv file with Clinvar labels
df = pd.read_csv("2023-03-06_figure4_gnomad_input.csv")
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

## set limits to the axis
ax.set_xlim([0, 0.05])
ax.set_ylim([0, 10000])

plt.tick_params(axis='both', which='major', labelsize=20)
fig.tight_layout()
plt.show()



