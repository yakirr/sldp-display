from __future__ import print_function, division
import numpy as np
import os
import pandas as pd
import pyutils.fs as fs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from plot import params, biaspower; reload(biaspower)

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/2.vary_h2v/compiled_results/'
outname = me+'/out/mainfig.biaspower.raw.pdf'

tickprops = {
        'direction':'out',
        'length':2,
        'width':0.8,
        'pad':4,
        'labelsize':6}
powererrorbarprops = {
        'linewidth':0.5,
        'markersize':1}
labelfontsize=6

desc='maf5'
refpanel='KG3.wim9nm'
estimands = {
        'mu':'$\mu$',
        'r_f':'$r_f$',
        'h2v':'$h^2_v$',
        'h2v_h2g':'$h^2_v/h^2_g$'
        }
weights = {
        'Winv_ahat_h':'non-trivial weights',
        'Winv_ahat_I':'trivial weights'
        }

## set up figure
fig = plt.figure(figsize=(5,2.5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])

## make bias plot for part a
print('making bias plot')
estimand='r_f'
biaspower.bias_plot(ax1, indir, desc, 'Winv_ahat_h', refpanel, estimand)
ax1.tick_params(**tickprops)
ax1.set_xlabel(r'True '+estimands[estimand], fontsize=labelfontsize)
ax1.set_ylabel(r'Estimated '+estimands[estimand], fontsize=labelfontsize)


## create power plot for part b
print('making power plot')
m='r_f'
hue, (r1,w1), (r2,w2) = ('Weight scheme',
            ('KG3.wim9nm','Winv_ahat_h'),('KG3.wim9nm','Winv_ahat_I'))

biaspower.power_plot(ax2,
        params.sldp+'/2.vary_h2v/compiled_results/',
        desc, 'Winv_ahat_h', refpanel,
        'r_f', '{:.2f}',
        label=weights['Winv_ahat_h'],
        c='b', labelfontsize=labelfontsize, **powererrorbarprops)
biaspower.power_plot(ax2,
        params.sldp+'/2.vary_h2v/compiled_results/',
        desc, 'Winv_ahat_I', refpanel,
        'r_f', '{:.2f}',
        label=weights['Winv_ahat_I'],
        c='r', labelfontsize=labelfontsize, **powererrorbarprops)
ax2.axis((-0.005, 0.055, 0, 1.05))
ax2.set_xlabel(r'True $r_f$', fontsize=labelfontsize)
ax2.tick_params(**tickprops)
ax2.legend(loc='upper left', fontsize=5, markerscale=2, borderpad=0.1,
        labelspacing=0.2, columnspacing=0.2)


# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print(m, ': saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
