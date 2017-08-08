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

# aesthetics
pie_props = {
        'startangle':180,
        'radius':1.4,
        'wedgeprops':{
            'linewidth':0.05,
            'edgecolor':'white',
            'width':0.4}
        }

# set parameters
nrows = 4; ncols = 10

# set up figure
fig = plt.figure(figsize=(6,3.5))
gs = gridspec.GridSpec(nrows,ncols)
ax1 = plt.subplot(gs[0,:7])
ax2 = plt.subplot(gs[1,:7])
ax4 = plt.subplot(gs[2,:7])
ax5 = plt.subplot(gs[3,:7])
ax3 = plt.subplot(gs[:2,7:])
ax6 = plt.subplot(gs[2:,7:])

# get data
global results, phenos
results = results_overview.init(
        [indir_bp+'/blueprint.maf5_Winv_ahat_h.all',
            indir_ntr+'/ntr.all'],
        [indir_bp+'/blueprint.maf5_Winv_ahat_h.fdr',
            indir_ntr+'/ntr.fdr'])
phenos = results[results.passed].pheno.unique()
summary = pd.concat([results_overview.summary_table(results, p) for p in phenos], axis=0)

passed = results[results.passed].copy()
passed['piecount'] = 1
passed['activator'] = passed.uniprot_activator# & ~passed.uniprot_repressor
passed['chromatin'] = False

# add info on validated chromatin modifiers in immune cells
valid_modifiers = [
        'EP300', # K27AC is main function
        'PU1', # K27AC, monocytes https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4705427/
        'PU1', # K4ME1, monocytes https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2898526/
                # http://www.sciencedirect.com/science/article/pii/S1097276510003667
        'CEBPB', # K4ME1, monocytes https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4039984/#R27
                # http://www.cell.com/stem-cell-reports/fulltext/S2213-6711(15)00188-5
                # http://www.cell.com/stem-cell-reports/fulltext/S2213-6711(16)30305-8
        # 'MYC', # causal arrow likely in other direction?
                # seems to create other histone marks, http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0003650
        # 'MAX', # causal arrow likely in other direction?
        'CEBPD', # K4ME1, monocytes
        'BATF' # T cells, http://www.nature.com/ni/journal/v18/n4/full/ni.3683.html?WT.feed_name=subjects_gene-regulation-in-immune-cells
        ]
for g in valid_modifiers:
    passed.loc[passed.genegroup==g, 'chromatin'] = True

myresults = passed[passed.pheno.str.contains('gene') & ~passed.pheno.str.contains('NTR')]
results_overview.segmented_bar(ax1, passed,
        ['BP_mono_gene_nor_combat_peer_10',
            'BP_neut_gene_nor_combat_peer_10'],
        'activator', 'r',
        'Expression (BLUEPRINT)', params.labelfontsize+2)

myresults = passed[passed.pheno.str.contains('NTR')]
results_overview.segmented_bar(ax2, passed,
        ['geneexp_total_NTR'],
        'activator', 'r',
        'Expression (NTR)', params.labelfontsize+2)

myresults = passed[passed.pheno.str.contains('K4ME1')]
results_overview.segmented_bar(ax4, myresults,
        ['BP_neut_K4ME1_log2rpm_peer_10',
            'BP_mono_K4ME1_log2rpm_peer_10'],
        'activator', 'r', #'chromatin', 'purple',
        'H3K4me1 (BLUEPRINT)', params.labelfontsize+2)

myresults = passed[passed.pheno.str.contains('K27AC')]
results_overview.segmented_bar(ax5, myresults,
        ['BP_neut_K27AC_log2rpm_peer_10'],
        'activator', 'r', #'chromatin', 'purple',
        'H3K27ac (BLUEPRINT)', params.labelfontsize+2)

# scatter plot for NTR vs BLUEPRINT
x = results[results.pheno.str.contains('neut') &
        results.pheno.str.contains('gene')].sort_values('annot')
y = results[results.pheno.str.contains('NTR')].sort_values('annot')
x['scatter_color'] = 'b'
x.loc[x.uniprot_activator & ~x.uniprot_repressor, 'scatter_color'] = 'r'
# x.loc[~x.uniprot_activator & x.uniprot_repressor, 'scatter_color'] = 'b'
mask = (x.scatter_color=='r').values
ax3.scatter(x[mask].sf_z, y[mask].sf_z, c=x[mask].scatter_color, alpha=0.7,
        s=2.5, linewidth=0, label='activators')
# mask = ((x.scatter_color=='b') & x.passed).values
# ax3.scatter(x[mask].sf_z, y[mask].sf_z, c=x[mask].scatter_color, alpha=1,
#         s=2.5, linewidth=0, label='repressors')
mask = (x.scatter_color!='r').values
ax3.scatter(x[mask].sf_z, y[mask].sf_z, c=x[mask].scatter_color, alpha=0.7,
        s=2.5, linewidth=0, label='other')
# mask = ~x.passed.values
# ax3.scatter(x[mask].sf_z, y[mask].sf_z, c=x[mask].scatter_color, alpha=0.5,
#         s=2.5, linewidth=0)
threshy = np.min(y[y.passed].sf_z)
threshx = np.min(x[x.passed].sf_z)
ax3.axhline(y=threshy, color='gray', linestyle='--', linewidth=0.5, alpha=0.8)
ax3.axvline(x=threshx, color='gray', linestyle='--', linewidth=0.5, alpha=0.8)
ax3.legend(fontsize=6, markerscale=2, borderpad=0.1,
    labelspacing=0.4, columnspacing=0.1, handletextpad=0, loc='upper left')
# add trendline
ax3.plot(np.unique(x.sf_z), np.poly1d(np.polyfit(x.sf_z, y.sf_z, 1))(np.unique(x.sf_z)),
        c='gray', linewidth=1.2)
r2 = np.corrcoef(x.sf_z, y.sf_z)[0,1]**2
ax3.text(3, -2, '$r^2$: {:.2g}'.format(r2),
        fontsize=params.labelfontsize, color='black')

ax3.set_xlabel(r'BLUEPRINT z-score', fontsize=params.labelfontsize)
ax3.set_ylabel(r'NTR z-score', fontsize=params.labelfontsize)
ax3.tick_params(**params.tickprops)

# add legend in bottom right
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
