from __future__ import print_function, division
import numpy as np
import os
import pandas as pd
import pyutils.fs as fs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from plot import params, biaspower

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/2.vary_h2v/compiled_results/'
outname = me+'/out/mainfig.biaspower.raw.pdf'
numerical_outname = me+'/out/xlsxtable.biaspower.xlsx'
writer = pd.ExcelWriter(numerical_outname)

# set aesthetics
powererrorbarprops = {
        'linewidth':0.5,
        'capsize':0,
        'markersize':1
        }

# set pretty text
estimands = {
        'mu':'$\mu$',
        'r_f':'$r_f$',
        'h2v':'$h^2_v$',
        'h2v_h2g':'$h^2_v/h^2_g$'
        }
weights = {
        'Winv_ahat_h':'generalized least-sq.',
        'Winv_ahat_I':'ordinary least-sq.',
        'RV':'naiveRV',
        'NA':'naive'
        }

# set parameters
desc='maf5'
refpanel='KG3.wim9nm'

## set up figure
fig = plt.figure(figsize=(5,2.5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])

## create power plot for part a
print('making power plot')
estimand='r_f'
numbers = biaspower.power_plot(ax1,
        params.sldp+'/2.vary_h2v/compiled_results/',
        desc, 'Winv_ahat_h', refpanel,
        estimand, '{:.3f}',
        label=weights['Winv_ahat_h'],
        c='b', labelfontsize=params.labelfontsize, **powererrorbarprops)
numbers.to_excel(writer, 'A1. power (SLDP)', index=False)

numbers = biaspower.power_plot(ax1,
        params.sldp+'/2.vary_h2v/compiled_results/',
        desc, 'Winv_ahat_I', refpanel,
        estimand, '{:.3f}',
        label=weights['Winv_ahat_I'],
        c='r', labelfontsize=params.labelfontsize, **powererrorbarprops)
numbers.to_excel(writer, 'A2. power (ordinary reg.)', index=False)
ax1.axis((-0.005, 0.055, 0, 1.05))
ax1.set_xlabel(r'True $r_f$', fontsize=params.labelfontsize)
ax1.tick_params(**params.tickprops)
ax1.legend(loc='upper left', fontsize=6, handlelength=0.8, borderpad=0.1,
        labelspacing=0.2, handletextpad=0.5)

## make bias plot for part b
print('making bias plot')
estimand='r_f'
numbers = biaspower.bias_plot(ax2, indir, desc, 'Winv_ahat_h', refpanel, estimand)
ax2.tick_params(**params.tickprops)
ax2.set_xlabel(r'True '+estimands[estimand], fontsize=params.labelfontsize)
ax2.set_ylabel(r'Estimated '+estimands[estimand], fontsize=params.labelfontsize)
numbers.to_excel(writer, 'B. bias and variance', index=False)



# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
writer.save()
