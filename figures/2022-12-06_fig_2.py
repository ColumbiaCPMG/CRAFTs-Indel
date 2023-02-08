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
gnomad_stats = pd.read_csv("2022-12-30_gnomad_fig2_input.csv")
percent_rare_sAF_rare_rAF = gnomad_stats['percent_rare_sAF_rare_rAF_gnomad']
percent_rare_sAF_common_rAF= gnomad_stats['percent_rare_sAF_common_rAF_gnomad']
rare_sAF_rare_rAF = gnomad_stats['rare_sAF_rare_rAF_gnomad']
rare_sAF_common_rAF = gnomad_stats['rare_sAF_common_rAF_gnomad']

## These are just manual inputs of the information in the csv file above 
#percent_rare_sAF_rare_rAF = ['86.8%', '81.9%', '78.2%', '75.2%']
#percent_rare_sAF_common_rAF= ['13.2%', '18.1%', '21.8%', '24.8%']
#rare_sAF_rare_rAF = [756541, 713070, 681345, 654913]
#rare_sAF_common_rAF = [114633, 158104, 189829, 216261]

"""gnomad section """
fig, ax = plt.subplots()
bar_x = labels
bar_height = rare_sAF_common_rAF
labels = ['10bps window','20bps window','30bps window','40bps window']
x = np.arange(len(labels))
width = 0.35
# bar_label = ['114633\n', '158104\n', '189829\n', '216261\n']
#bar_label = ['13%', '18%', '22%', '25%']
bar_label_rounded = [round(num) for num in percent_rare_sAF_common_rAF.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot = plt.bar(x - width/2  ,bar_height,width = 0.35,
                   label = 'gnomAD: \nrare sAF and common rAF',color = [(0, 0, 0)])
autolabel_bot(bar_plot)
# bar_label = ['756,541\n', '713,070\n', '681,345\n', '654,913\n']
#bar_label = ['87%', '82%', '78%', '75%']
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
IGM_stats = pd.read_csv("2022-12-30_IGM_fig2_input.csv")
percent_rare_sAF_rare_rAF = IGM_stats['percent_rare_sAF_rare_rAF_IGM']
percent_rare_sAF_common_rAF= IGM_stats['percent_rare_sAF_common_rAF_IGM']
rare_sAF_rare_rAF = IGM_stats['rare_sAF_rare_rAF_IGM']
rare_sAF_common_rAF = IGM_stats['rare_sAF_common_rAF_IGM']

## These are just manual inputs of the information in the csv file above 
#percent_rare_sAF_common_rAF = ['16.1%', '22.7%', '27.3%', '30.8%']
#percent_rare_sAF_rare_rAF = ['83.9%', '77.3%', '72.7%', '69.2%']
#rare_sAF_common_rAF = [40496, 57157, 68719, 77528]
#rare_sAF_rare_rAF = [211142, 194481, 182919, 174110]

'''IGM section'''
bar_height = rare_sAF_common_rAF
labels = ['10bps window','20bps window','30bps window','40bps window']
x = np.arange(len(labels))
witdth = 0.35
# bar_label = ['40,496\n', '57,157\n', '68,719\n', '77,528\n']
#bar_label = ['16%', '23%', '27%', '31%']
bar_label_rounded = [round(num) for num in percent_rare_sAF_common_rAF.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot = plt.bar(x + width/2 ,bar_height,width = 0.35,
                   label = 'IGM: \nrare sAF and common rAF',color = [(0.25,0.25,0.25)])
autolabel_bot(bar_plot)

# bar_label = ['211,142\n', '194,481\n', '182,919\n', '174,110\n']
#bar_label = ['84%', '77%', '73%', '69%']
bar_label_rounded = [round(num) for num in percent_rare_sAF_rare_rAF.tolist()]
bar_label = [str(num) + "%" for num in bar_label_rounded]
bar_plot2 = plt.bar(x + width/2,rare_sAF_rare_rAF,width = 0.35,bottom = rare_sAF_common_rAF,color = [0.75,0.75,0.75],
                    label = 'IGM \nrare sAF and rare rAF')
autolabel_top(bar_plot2)
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()
