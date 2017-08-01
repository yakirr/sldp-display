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
outname = me+'/out/suppfig.power_refpanel.pdf'

# set aesthetics
powererrorbarprops = {
        'linewidth':1}

# set pretty text
refpanels = {
        'KG3.wim9nm' : '1KG reference ($N=489$)',
        'GERAimp.wim9nm.ref500' : 'GERA reference ($N=500$)'
        }

# set parameters
desc='maf5'
weights='Winv_ahat_h'

# set up figure
fig = plt.figure(figsize=(2.5,2.5))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])

## create power plot for part b
print('making power plot')
estimand = 'r_f'
weights = 'Winv_ahat_h'
biaspower.power_plot(ax1,
        params.sldp+'/2.vary_h2v/compiled_results/',
        desc, weights, 'KG3.wim9nm',
        estimand, '{:.2f}',
        label=refpanels['KG3.wim9nm'],
        c='b', labelfontsize=params.labelfontsize, **powererrorbarprops)
biaspower.power_plot(ax1,
        params.sldp+'/2.vary_h2v/compiled_results/',
        desc, weights, 'GERAimp.wim9nm.ref500',
        estimand, '{:.2f}',
        label=refpanels['GERAimp.wim9nm.ref500'],
        c='g', labelfontsize=params.labelfontsize, **powererrorbarprops)
ax1.axis((-0.005, 0.055, 0, 1.05))
ax1.set_xlabel(r'True $r_f$', fontsize=params.labelfontsize)
ax1.tick_params(**params.tickprops)
ax1.legend(loc='upper left', fontsize=5, markerscale=2, borderpad=0.1,
        labelspacing=0.2, columnspacing=0.2)


# finishing touches and save
sns.despine()
plt.tight_layout()

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
