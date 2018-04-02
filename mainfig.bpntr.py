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
indir = params.sldprev+'/1.basset1tfs_p12/'
# outname = me+'/out/mainfig.bpntr.raw.pdf'
outname = '/n/scratch2/yar2/mainfig.bpntr.raw.pdf'

## aesthetics
scatterprops_r = {
        'alpha':0.8,
        's':2.5,
        'linewidth':0
        }
scatterprops_b = {
        'alpha':0.5,
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
        [indir+'/molecular_NTR/all.gwresults', indir+'/molecular_BP/all.gwresults'],
        [indir+'/molecular_NTR/fdr5.gwresults', indir+'/molecular_BP/fdr5.gwresults'],
        'desc')
results['activator'] = results.uniprot_activator & ~results.uniprot_repressor
results['repressor'] = ~results.uniprot_activator & results.uniprot_repressor
results['ambig'] = ~results.activator & ~results.repressor
results['activatingness'] = 2*results.activator.astype(int) + results.ambig.astype(int)
passed = results[results.passed].copy()

## make figure
# BP gene exp
myresults = passed[passed.pheno.str.contains('gene') & ~passed.pheno.str.contains('NTR')]
print(myresults.pheno.value_counts())
nassoc = len(myresults)
myresults = myresults.sort_values(by='p_fast').iloc[:100]
results_overview.segmented_bar(ax1, myresults,
        ['BPh2g50sqrt_mono_gene_nor_combat',
            'BPh2g50sqrt_neut_gene_nor_combat',
            'BPh2g50sqrt_tcel_gene_nor_combat'],
        {'activator':'r', 'repressor':'b', 'ambig':'#aaaaaa'},
        'Expression, BLUEPRINT (top 100 of {} associations)'.format(nassoc),
        params.labelfontsize)

# NTR gene exp
myresults = passed[passed.pheno.str.contains('NTR')]
print(myresults.pheno.value_counts())
nassoc = len(myresults)
myresults = myresults.sort_values(by='p_fast').iloc[:100]
results_overview.segmented_bar(ax2, myresults,
        ['NTRh2g50sqrt_blood_gene'],
        {'activator':'r', 'repressor':'b', 'ambig':'#aaaaaa'},
        'Expression, NTR ({} associations)'.format(nassoc),
        params.labelfontsize)

# BP K4me1
myresults = passed[passed.pheno.str.contains('K4ME1')]
print(myresults.pheno.value_counts())
nassoc = len(myresults)
myresults = myresults.sort_values(by='p_fast').iloc[:100]
results_overview.segmented_bar(ax4, myresults,
        ['BPh2g50sqrt_neut_K4ME1_log2rpm',
            'BPh2g50sqrt_mono_K4ME1_log2rpm',
            'BPh2g50sqrt_tcel_K4ME1_log2rpm'],
        {'activator':'r', 'repressor':'b', 'ambig':'#aaaaaa'},
        'H3K4me1, BLUEPRINT (top 100 of {} associations)'.format(nassoc),
        params.labelfontsize)

# BP K27ac
myresults = passed[passed.pheno.str.contains('K27AC')]
print(myresults.pheno.value_counts())
nassoc = len(myresults)
myresults = myresults.sort_values(by='p_fast').iloc[:100]
results_overview.segmented_bar(ax5, myresults,
        ['BPh2g50sqrt_neut_K27AC_log2rpm',
            'BPh2g50sqrt_mono_K27AC_log2rpm',
            'BPh2g50sqrt_tcel_K27AC_log2rpm'],
        {'activator':'r', 'repressor':'b', 'ambig':'#aaaaaa'},
        'H3K27ac, BLUEPRINT (top 100 of {} associations)'.format(nassoc),
        params.labelfontsize)

## add global figure legend in bottom right
patches = results_overview.legend_contents('desc')
ax6.set_axis_off()
ax6.legend(handles=patches, fontsize=7, handlelength=1.5, borderpad=0.1,
    labelspacing=0.4, loc='upper left', bbox_to_anchor=(-0.18,1))

## scatter plot for NTR vs BLUEPRINT
# create data
x = results[results.pheno.str.contains('neut') &
        results.pheno.str.contains('gene')].sort_values('annot')
y = results[results.pheno.str.contains('NTR')].sort_values('annot')

# plot points
ax3.scatter(x[x.activator.values].z_fast, y[x.activator.values].z_fast, c='r',
        label='activating',
        **scatterprops_r)
ax3.scatter(x[x.ambig.values].z_fast, y[x.ambig.values].z_fast, c='#aaaaaa',
        label='ambiguous',
        **scatterprops_b)
ax3.scatter(x[x.repressor.values].z_fast, y[x.repressor.values].z_fast, c='b',
        label='repressing',
        **scatterprops_b)

# add horizontal and vertical significance threshold lines
threshy = np.min(y[y.passed].z_fast)
threshx = np.min(x[x.passed].z_fast)
ax3.axhline(y=threshy, **params.sig_thresh_line_props)
ax3.axvline(x=threshx, **params.sig_thresh_line_props)

# add trendline
ax3.plot(np.unique(x.z_fast), np.poly1d(np.polyfit(x.z_fast, y.z_fast, 1))(np.unique(x.z_fast)),
        **trendlineprops)

# add r2 text label
r = np.corrcoef(x.z_fast, y.z_fast)[0,1]
ax3.text(3.5, -1.5, r'$r = {:.2g}$'.format(r),
        fontsize=params.labelfontsize, color='black')

# add legend and format axes
ax3.legend(fontsize=6, markerscale=1.5, borderpad=0.05,
    labelspacing=0.1, columnspacing=0.01, handletextpad=-0.5, loc=(-0.03,0.75))
ax3.set_xlim(-3, 7)
ax3.set_ylim(-2.1, 5)
ax3.set_xticks([-2,0,2,4,6])
ax3.set_yticks([-2,0,2,4])
ax3.set_xlabel(r'BLUEPRINT z-score', fontsize=params.labelfontsize)
ax3.set_ylabel(r'NTR z-score', fontsize=params.labelfontsize)
ax3.tick_params(**params.tickprops)

# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
