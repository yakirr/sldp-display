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
outname = me+'/out/suppfig.bias_moreestimands.pdf'

# set aesthetics

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
weights='Winv_ahat_h'

## set up figure
fig = plt.figure(figsize=(5,5))
gs = gridspec.GridSpec(2,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[1,0])
ax4 = plt.subplot(gs[1,1])

# make bias plots for all four estimands
print('making bias plot')
for ax, (e, name) in zip([ax1,ax2,ax3,ax4], estimands.items()):
    print(e, name)
    biaspower.bias_plot(ax, indir, desc, weights, refpanel, e)
    ax.tick_params(**dict(params.tickprops, labelsize=6))
    ax.set_xlabel(r'True '+name, fontsize=params.labelfontsize)
    ax.set_ylabel(r'Estimated '+name, fontsize=params.labelfontsize)


# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
