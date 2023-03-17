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
        #ax.text(rect.get_x() + rect.get_width()/2., 0.35*height,
        ax.text(rect.get_x() + rect.get_width()/2., 0.25*height,
                bar_label[idx],
                ha='center', va='bottom', rotation=360, color = 'white')

def autolabel_top(rects):
    for idx,rect in enumerate(bar_plot):
        height = rect.get_height()
        #ax.text(rect.get_x() + rect.get_width()/2., 3.5*height,
        ax.text(rect.get_x() + rect.get_width()/2., 2.5 *height,
                bar_label[idx],
                ha='center', va='bottom', rotation=360, color = 'white')

labels = ['10bps window','20bps window','30bps window','40bps window']

'''gnomad stats'''
gnomad_stats = pd.read_csv("2023-03-06_gnomad_fig2_input.csv")
percent_rare_sAF_rare_rAF = gnomad_stats['percent_rare_sAF_rare_rAF_gnomad']
percent_rare_sAF_common_rAF= gnomad_stats['percent_rare_sAF_common_rAF_gnomad']
rare_sAF_rare_rAF = gnomad_stats['rare_sAF_rare_rAF_gnomad']
rare_sAF_common_rAF = gnomad_stats['rare_sAF_common_rAF_gnomad']

"""gnomad section """
fig, ax = plt.subplots()
bar_x = labels
bar_height = rare_sAF_common_rAF
labels = ['10bps window','20bps window','30bps window','40bps window']
x = np.arange(len(labels))
width = 0.35
bar_label_rounded = [round(num) for num in percent_rare_sAF_common_rAF.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot = plt.bar(x - width/2  ,bar_height,width = 0.35,
                   label = 'gnomAD: \nrare sAF and common rAF',color = [(0, 0, 0)])
autolabel_bot(bar_plot)
bar_label_rounded = [round(num) for num in percent_rare_sAF_rare_rAF.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot2 = plt.bar(x - width/2,rare_sAF_rare_rAF,width = 0.35,bottom = rare_sAF_common_rAF,
                    color = [(0.5,0.5,0.5)],label = 'gnomAD: \nrare sAF and rare rAF')
autolabel_top(bar_plot2)
ax.set_xticks(x, labels)
plt.title('Indels that are rare by sAF and common by rAF/ indels that are rare by both sAF and rAF \n for IGM and gnomAD dataset')
plt.ylabel('Number of Indels')
plt.xlabel('Sliding Windows')
plt.legend()

'''IGM stats'''
IGM_stats = pd.read_csv("2023-03-06_IGM_fig2_input.csv")
percent_rare_sAF_rare_rAF = IGM_stats['percent_rare_sAF_rare_rAF_IGM']
percent_rare_sAF_common_rAF= IGM_stats['percent_rare_sAF_common_rAF_IGM']
rare_sAF_rare_rAF = IGM_stats['rare_sAF_rare_rAF_IGM']
rare_sAF_common_rAF = IGM_stats['rare_sAF_common_rAF_IGM']


'''IGM section'''
bar_height = rare_sAF_common_rAF
labels = ['10bps window','20bps window','30bps window','40bps window']
x = np.arange(len(labels))
width = 0.35

bar_label_rounded = [round(num) for num in percent_rare_sAF_common_rAF.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot = plt.bar(x + width/2 ,bar_height,width = 0.35,
                   label = 'IGM: \nrare sAF and common rAF',color = [(0.25,0.25,0.25)])
autolabel_bot(bar_plot)


bar_label_rounded = [round(num) for num in percent_rare_sAF_rare_rAF.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot2 = plt.bar(x + width/2,rare_sAF_rare_rAF,width = 0.35,bottom = rare_sAF_common_rAF,color = [0.75,0.75,0.75],
                    label = 'IGM \nrare sAF and rare rAF')
autolabel_top(bar_plot2)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
