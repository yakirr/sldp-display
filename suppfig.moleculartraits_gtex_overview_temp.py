from __future__ import print_function, division
import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.cluster.hierarchy as sch
import ypy.fs as fs
from plot import params, results_overview; reload(results_overview)

me = os.path.dirname(os.path.abspath(__file__))
gtexindir = params.sldp+'../newsldp/7.process_gtex/compiled_results/'
oldindir = params.sldp+'/7.p9_a9/compiled_results/'
# outname = me+'/out/mainfig.moleculartraits_overview.raw.pdf'

## aesthetics
scatterprops = {
        'alpha':0.7,
        's':2.5,
        'linewidth':0
        }
trendlineprops = {
        'c':'gray',
        'linewidth':1.2
        }

## set up figure
fig = plt.figure(figsize=(6,3.5))
gs = gridspec.GridSpec(4,10)
ax1 = plt.subplot(gs[0,:7])
ax2 = plt.subplot(gs[1,:7])
ax4 = plt.subplot(gs[2,:7])
ax5 = plt.subplot(gs[3,:7])
ax3 = plt.subplot(gs[:2,7:])
ax6 = plt.subplot(gs[2:,7:])

## get data
global results, phenos
results = results_overview.init(
        [gtexindir+'/gtex.gwresults', oldindir+'/p9.molecular.all'],
        [gtexindir+'/gtex.fdr5', oldindir+'/p9.molecular.fdr5'],
        'desc')
results['activator'] = results.uniprot_activator & ~results.uniprot_repressor
passed = results[results.passed].copy()

## make figure
## scatter plot for comparison of two data sets
# create data
x = results[results.pheno.str.contains('gtex') &
        results.pheno.str.contains('Whole_Blood')].sort_values('annot')
y = results[results.pheno.str.contains('gtex') &
        results.pheno.str.contains('Nerve_Tibial')].sort_values('annot')
# y = y.drop(['z'], axis=1).rename(columns={'sf_z':'z', 'sf_p':'p'}) for NTR or BLUEPRINT

# plot points
ax3.scatter(x[x.activator.values].z, y[x.activator.values].z, c='r',
        label='activating',
        **scatterprops)
ax3.scatter(x[~x.activator.values].z, y[~x.activator.values].z, c='b',
        label='other',
        **scatterprops)

# add horizontal and vertical significance threshold lines
threshy = np.min(y[y.passed].z)
threshx = np.min(x[x.passed].z)
ax3.axhline(y=threshy, **params.sig_thresh_line_props)
ax3.axvline(x=threshx, **params.sig_thresh_line_props)

# add trendline
ax3.plot(np.unique(x.z), np.poly1d(np.polyfit(x.z, y.z, 1))(np.unique(x.z)),
        **trendlineprops)

# add r2 text label
r = np.corrcoef(x.z, y.z)[0,1]
ax3.text(3, -2, '$r = {:.2g}$'.format(r),
        fontsize=params.labelfontsize, color='black')

# add legend and format axes
ax3.legend(fontsize=6, markerscale=2, borderpad=0.1,
    labelspacing=0.4, columnspacing=0.1, handletextpad=0, loc='upper left')
ax3.set_xticks([-2,0,2,4])
ax3.set_yticks([-2,0,2,4])
ax3.set_xlabel(r'GTEx z-score', fontsize=params.labelfontsize)
ax3.set_ylabel(r'NTR z-score', fontsize=params.labelfontsize)
ax3.tick_params(**params.tickprops)

# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print('saving figure')
plt.show()
# fs.makedir_for_file(outname)
# plt.savefig(outname)
# plt.close()
