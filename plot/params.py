from __future__ import print_function, division
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='ticks', palette='Set2')
plt.rc('axes', linewidth=0.8)

sldp = '/groups/price/yakir/sldp/'

# aesthetics
labelfontsize=6
tickprops = {
        'direction':'out',
        'length':2,
        'width':0.8,
        'pad':4,
        'labelsize':6}
sig_thresh_line_props = {
        'color':'gray',
        'linestyle':'--',
        'linewidth':0.5,
        'alpha':0.8
        }
qqprops = {
        's':1.5,
        'fontsize':labelfontsize,
        'linewidth':0}
