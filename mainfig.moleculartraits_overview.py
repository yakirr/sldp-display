from __future__ import print_function, division
import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.cluster.hierarchy as sch
import pyutils.fs as fs
from plot import params, results_overview; reload(results_overview)

me = os.path.dirname(os.path.abspath(__file__))
indir_bp = params.sldp+'/7.blueprint_a9/compiled_results_uniform_normalization/'
indir_ntr = params.sldp+'/8.totalexpNTR_a9/compiled_results_uniform_normalization/'
outname = me+'/out/mainfig.moleculartraits_overview.raw.pdf'

## aesthetics
scatterprops = {
        'alpha':0.7,
        's':2.5,
        'linewidth':0
        }
sig_thresh_line_props = {
        'color':'gray',
        'linestyle':'--',
        'linewidth':0.5,
        'alpha':0.8
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
        [indir_bp+'/blueprint.maf5_Winv_ahat_h.all',
            indir_ntr+'/ntr.all'],
        [indir_bp+'/blueprint.maf5_Winv_ahat_h.fdr',
            indir_ntr+'/ntr.fdr'])
passed = results[results.passed].copy()
passed['activator'] = passed.uniprot_activator# & ~passed.uniprot_repressor

## make figure
# BP gene exp
myresults = passed[passed.pheno.str.contains('gene') & ~passed.pheno.str.contains('NTR')]
results_overview.segmented_bar(ax1, passed,
        ['BP_mono_gene_nor_combat_peer_10',
            'BP_neut_gene_nor_combat_peer_10'],
        'activator', 'r',
        'Expression (BLUEPRINT)', params.labelfontsize+2)

# NTR gene exp
myresults = passed[passed.pheno.str.contains('NTR')]
results_overview.segmented_bar(ax2, passed,
        ['geneexp_total_NTR'],
        'activator', 'r',
        'Expression (NTR)', params.labelfontsize+2)

# BP K4me1
myresults = passed[passed.pheno.str.contains('K4ME1')]
results_overview.segmented_bar(ax4, myresults,
        ['BP_neut_K4ME1_log2rpm_peer_10',
            'BP_mono_K4ME1_log2rpm_peer_10'],
        'activator', 'r',
        'H3K4me1 (BLUEPRINT)', params.labelfontsize+2)

# BP K27ac
myresults = passed[passed.pheno.str.contains('K27AC')]
results_overview.segmented_bar(ax5, myresults,
        ['BP_neut_K27AC_log2rpm_peer_10'],
        'activator', 'r',
        'H3K27ac (BLUEPRINT)', params.labelfontsize+2)

## scatter plot for NTR vs BLUEPRINT
# create data
x = results[results.pheno.str.contains('neut') &
        results.pheno.str.contains('gene')].sort_values('annot')
y = results[results.pheno.str.contains('NTR')].sort_values('annot')

# plot points
ax3.scatter(x[x.activator].sf_z, y[x.activator].sf_z, c='r',
        label='activators',
        **scatterprops)
ax3.scatter(x[~x.activator].sf_z, y[~x.activator].sf_z, c='b',
        label='other',
        **scatterprops)

# add horizontal and vertical significance threshold lines
threshy = np.min(y[y.passed].sf_z)
threshx = np.min(x[x.passed].sf_z)
ax3.axhline(y=threshy, **sig_thresh_line_props)
ax3.axvline(x=threshx, **sig_thresh_line_props)

# add trendline
ax3.plot(np.unique(x.sf_z), np.poly1d(np.polyfit(x.sf_z, y.sf_z, 1))(np.unique(x.sf_z)),
        **trendlineprops)

# add r2 text label
r2 = np.corrcoef(x.sf_z, y.sf_z)[0,1]**2
ax3.text(3, -2, '$r^2$: {:.2g}'.format(r2),
        fontsize=params.labelfontsize, color='black')

# add legend and format axes
ax3.legend(fontsize=6, markerscale=2, borderpad=0.1,
    labelspacing=0.4, columnspacing=0.1, handletextpad=0, loc='upper left')
ax3.set_xticks([-2,0,2,4])
ax3.set_yticks([-2,0,2,4])
ax3.set_xlabel(r'BLUEPRINT z-score', fontsize=params.labelfontsize)
ax3.set_ylabel(r'NTR z-score', fontsize=params.labelfontsize)
ax3.tick_params(**params.tickprops)

## add global figure legend in bottom right
patches = results_overview.legend_contents('desc')
ax6.set_axis_off()
ax6.legend(handles=patches, fontsize=6, markerscale=2, borderpad=0.1,
    labelspacing=0.4, columnspacing=0.2, loc='upper left')

# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
