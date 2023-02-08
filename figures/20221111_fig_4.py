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
gnomad_stats = pd.read_csv("2023-01-06_gnomad_fig4_input.csv")
percent_in_LCR = gnomad_stats["percent_in_LCR"]
percent_outside_LCR = gnomad_stats["percent_outside_LCR"]
in_LCR = gnomad_stats["in_LCR"]
outside_LCR = gnomad_stats["outside_LCR"]

## These are manual inputs of the same info. in the csv file 
#percent_in_LCR= ['24.7%', '20.6%', '18.2%', '16.2%']
#percent_outside_LCR = ['75.3%', '79.4%', '81.8%', '83.8%']
#in_LCR = [28314, 32588, 34518, 35112]
#outside_LCR = [86319, 125516, 155311, 181149]

"""gnomad section """
fig, ax = plt.subplots()
bar_x = labels
bar_height = in_LCR
labels = ['10bps window','20bps window','30bps window','40bps window']
x = np.arange(len(labels))
width = 0.35
#bar_label = ['25%', '21%', '18%', '16%']
bar_label_rounded = [round(num) for num in percent_in_LCR.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot = plt.bar(x - width/2  ,bar_height,width = 0.35,
                   label = 'gnomAD: in LCR',color = [(0, 0, 0)])
autolabel_bot(bar_plot)
#bar_label = ['75%', '79%', '82%', '84%']
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
IGM_stats = pd.read_csv("2023-01-06_IGM_fig4_input.csv")
percent_in_LCR = IGM_stats["percent_in_LCR"]
percent_outside_LCR = IGM_stats["percent_outside_LCR"]
in_LCR = IGM_stats["in_LCR"]
outside_LCR = IGM_stats["outside_LCR"]

## These are manual inputs of the same info. in the csv file 
#percent_in_LCR = ['19.3%', '14.6%', '12.5%', '11.2%']
#percent_not_in_LCR = ['80.7%', '85.4%', '87.49%', '88.8%']
#in_LCR = [7799, 8358, 8595, 8695]
#outside_LCR = [32697,48799,60124,68833]

'''IGM section'''
bar_height = in_LCR
labels = ['10bps window','20bps window','30bps window','40bps window']
x = np.arange(len(labels))
width = 0.35
#bar_label = ['19%', '15%', '13%', '11%']
bar_label_rounded = [round(num) for num in percent_in_LCR.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot = plt.bar(x + width/2 ,bar_height,width = 0.35,
                   label = 'IGM: in LCR',color = [(0.25,0.25,0.25)])
autolabel_bot(bar_plot)

#bar_label = ['81%', '85%', '87%', '89%']
bar_label_rounded = [round(num) for num in percent_outside_LCR.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot2 = plt.bar(x + width/2,outside_LCR,width = 0.35,bottom = in_LCR,color = [0.75,0.75,0.75],
                    label = 'IGM: not in LCR')
autolabel_top(bar_plot2)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#plt.savefig('plot.jpg')
plt.show()
