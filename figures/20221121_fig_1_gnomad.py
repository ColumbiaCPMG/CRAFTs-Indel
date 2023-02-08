import seaborn as sn
import matplotlib
import matplotlib.pyplot as plt
import pickle
a3 = pickle.load(open('2022-12-07_gnomad_fig1_input.pkl', 'rb'))
bps = [5,10,15,20]
for bp in bps:
    plt.figure()
    sn.distplot(a3[bp],norm_hist = False, color = "gray")
    plt.xlabel('Bin Length (bps)', fontsize=15)
    #plt.ylabel('number of bins')
    plt.ylabel('Density', fontsize=15)
    plt.title('Bin Length Distribution for bp Window =' + str(bp * 2), fontsize=15)
    plt.tick_params(axis='both', which='major', labelsize=15)
    plt.xlim([0, 1150])
    plt.show()
