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
outname = me+'/out/suppfig.power_naive.pdf'

# set aesthetics
powererrorbarprops = {
        'linewidth':0.5,
        'capsize':0,
        'markersize':1}

# set pretty text
estimands = {
        'mu':'$\mu$',
        'r_f':'$r_f$',
        'h2v':'$h^2_v$',
        'h2v_h2g':'$h^2_v/h^2_g$'
        }

# set parameters
desc='maf5'
refpanel='KG3.wim9nm'
estimand='r_f'

## set up figure
fig = plt.figure(figsize=(5,2.5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])

## create power plot for part a
print('making power plot a')
biaspower.power_plot(ax1,
        params.sldp+'/2.vary_h2v/compiled_results/',
        'maf5', 'Winv_ahat_h', refpanel,
        estimand, '{:.3f}',
        label='SLDP (GLS)',
        c='b', labelfontsize=params.labelfontsize, **powererrorbarprops)
biaspower.power_plot(ax1,
        params.sldp+'/2.vary_h2v/compiled_results/',
        'maf5', 'Winv_ahat_I', refpanel,
        estimand, '{:.3f}',
        label='SLDP (OLS)',
        c='r', labelfontsize=params.labelfontsize, **powererrorbarprops)
biaspower.power_plot(ax1,
        params.sldp+'/2.vary_h2v/compiled_results/',
        'naive', 'Winv_ahat_I', refpanel,
        estimand, '{:.3f}',
        label='naive method',
        c='g', labelfontsize=params.labelfontsize, **powererrorbarprops)
ax1.axis((-0.005, 0.055, 0, 1.05))
ax1.set_xlabel(r'True $r_f$', fontsize=params.labelfontsize)
ax1.tick_params(**params.tickprops)
ax1.legend(loc='upper left', fontsize=5, markerscale=2, borderpad=0.1,
        labelspacing=0.2, columnspacing=0.2)

## create power plot for part b
print('making power plot b')
suff = ['_halfhm3_'+str(i) for i in range(1,6)]
biaspower.power_plot(ax2,
        params.sldp+'/2.vary_h2v/compiled_results/',
        'maf5',
        ['Winv_ahat_h'+s for s in suff],
        refpanel,
        estimand, '{:.3f}',
        samplesizefactor=5,
        label='SLDP (GLS)',
        c='b', labelfontsize=params.labelfontsize, **powererrorbarprops)
biaspower.power_plot(ax2,
        params.sldp+'/2.vary_h2v/compiled_results/',
        'maf5',
        ['Winv_ahat_I'+s for s in suff],
        refpanel,
        estimand, '{:.3f}',
        samplesizefactor=5,
        label='SLDP (OLS)',
        c='r', labelfontsize=params.labelfontsize, **powererrorbarprops)
biaspower.power_plot(ax2,
        params.sldp+'/2.vary_h2v/compiled_results/',
        'naive',
        ['Winv_ahat_I'+s for s in suff],
        refpanel,
        estimand, '{:.3f}',
        samplesizefactor=5,
        label='naive method',
        c='g', labelfontsize=params.labelfontsize, **powererrorbarprops)
ax2.axis((-0.005, 0.055, 0, 1.05))
ax2.set_xlabel(r'True $r_f$', fontsize=params.labelfontsize)
ax2.tick_params(**params.tickprops)
ax2.legend(loc='upper left', fontsize=5, markerscale=2, borderpad=0.1,
        labelspacing=0.2, columnspacing=0.2)

# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
