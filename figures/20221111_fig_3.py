## This code is to generate Figure 4 - LCR
'''scripts to make side by side stacked bar plots with number and percentage on each bar'''
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}
matplotlib.rc('font', **font)
def autolabel_bot(rects):
    for idx,rect in enumerate(bar_plot):
        height = rect.get_height()
        print(height)
        ax.text(rect.get_x() + rect.get_width()/2., 0.1*height,
                bar_label[idx],
                ha='center', va='bottom', rotation=360, color = 'white')

def autolabel_top(rects):
    for idx,rect in enumerate(bar_plot):
        height = rect.get_height()
        #ax.text(rect.get_x() + rect.get_width()/2., 3.5*height,
        ax.text(rect.get_x() + rect.get_width()/2., 1.1*height,
                bar_label[idx],
                ha='center', va='bottom', rotation=360, color = 'white')

labels = ['10bps window','20bps window','30bps window','40bps window']

'''gnomad stats'''
gnomad_stats = pd.read_csv("2023-03-17_gnomad_fig3_input.csv")
percent_in_LCR = gnomad_stats["percent_in_LCR"]
percent_outside_LCR = gnomad_stats["percent_outside_LCR"]
in_LCR = gnomad_stats["in_LCR"]
outside_LCR = gnomad_stats["outside_LCR"]

"""gnomad section """
fig, ax = plt.subplots()
bar_x = labels
bar_height = in_LCR
labels = ['10bps window','20bps window','30bps window','40bps window']
x = np.arange(len(labels))
width = 0.35
bar_label_rounded = [round(num) for num in percent_in_LCR.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot = plt.bar(x - width/2  ,bar_height,width = 0.35,
                   label = 'gnomAD: in LCR',color = [(0, 0, 0)])
autolabel_bot(bar_plot)
bar_label_rounded = [round(num) for num in percent_outside_LCR.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot2 = plt.bar(x - width/2,outside_LCR,width = 0.35,bottom = in_LCR,
                    color = [(0.5,0.5,0.5)],label = 'gnomAD: not in LCR')
autolabel_top(bar_plot2)
ax.set_xticks(x, labels)
plt.title('Distribution of rare sAF and common rAF indels inside and outside of low complexity regions (LCRs) in the gnomAD and IGM data sets.')
plt.ylabel('Number of Indels')
plt.xlabel('Sliding Windows')
plt.legend()

'''IGM stats'''
IGM_stats = pd.read_csv("2023-03-17_IGM_fig3_input.csv")
percent_in_LCR = IGM_stats["percent_in_LCR"]
percent_outside_LCR = IGM_stats["percent_outside_LCR"]
in_LCR = IGM_stats["in_LCR"]
outside_LCR = IGM_stats["outside_LCR"]


'''IGM section'''
bar_height = in_LCR
labels = ['10bps window','20bps window','30bps window','40bps window']
x = np.arange(len(labels))
width = 0.35
bar_label_rounded = [round(num) for num in percent_in_LCR.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot = plt.bar(x + width/2 ,bar_height,width = 0.35,
                   label = 'IGM: in LCR',color = [(0.25,0.25,0.25)],hatch="///")
autolabel_bot(bar_plot)

bar_label_rounded = [round(num) for num in percent_outside_LCR.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot2 = plt.bar(x + width/2,outside_LCR,width = 0.35,bottom = in_LCR,color = [0.75,0.75,0.75], hatch="///",
                    label = 'IGM: not in LCR')
autolabel_top(bar_plot2)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
