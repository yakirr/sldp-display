from __future__ import print_function, division
import os
import numpy as np
import pandas as pd
import pyutils.fs as fs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from plot import params, biaspower

me = os.path.dirname(os.path.abspath(__file__))
outname = me+'/out/suppfig.power_N_h2g.raw.pdf'

tickprops = {
        'direction':'out',
        'length':2,
        'width':0.8,
        'pad':4,
        'labelsize':6}
powererrorbarprops = {
        'c':'b',
        'linewidth':1}
labelfontsize=6

desc='maf5'
weights='Winv_ahat_h'
refpanel='KG3.wim9nm'

# set up figure
fig = plt.figure(figsize=(5,2.5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])

## power as a function of N
biaspower.power_plot(ax1,
        params.sldp+'/4.vary_N/compiled_results/',
        desc, weights, refpanel,
        'N', '{:.0f}',
        labelfontsize=labelfontsize,
        **powererrorbarprops)
ax1.axis((0, 50000, 0, 1))
ax1.set_xlabel('Sample size', fontsize=7)
ax1.tick_params(**tickprops)

## power as a function of h2g
biaspower.power_plot(ax2,
        params.sldp+'/3.vary_h2g/compiled_results/',
        desc, weights, refpanel,
        'h2g', '{:.2f}',
        labelfontsize=labelfontsize,
        **powererrorbarprops)
ax2.axis((0, 0.6, 0, 1))
ax2.set_xlabel('$h^2_g$', fontsize=7)
ax2.tick_params(**tickprops)


# finishing touches and save
sns.despine()
plt.tight_layout()

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
