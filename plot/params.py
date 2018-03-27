from __future__ import print_function, division
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', palette='Set2')
plt.rc('axes', linewidth=0.8)

sldp = '../sldp.sub_ng/'
sldprev = '../sldp.rev1_ng/'
sumstats = '../data/sumstats.hm3/'
kg3 = '../../ldsc/reference_files/1000G_EUR_Phase3/'

# aesthetics
labelfontsize=8
tickprops = {
        'direction':'out',
        'length':2,
        'width':0.8,
        'pad':4,
        'labelsize':7}
sig_thresh_line_props = {
        'color':'gray',
        'linestyle':'--',
        'linewidth':0.5,
        'alpha':0.8
        }
enrichment_thresh_line_props = {
        'color':'darkgray',
        'linestyle':'-',
        'linewidth':1,
        'alpha':1
        }
qqprops = {
        's':1.5,
        'fontsize':labelfontsize,
        'linewidth':0}
