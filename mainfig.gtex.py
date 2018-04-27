from __future__ import print_function, division
import pandas as pd
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats as st
import ypy.fs as fs
from plot import params, info

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldprev+'/1.basset1tfs_p12/old_results/molecular_gtexv7_tissues/'
outname = me+'/out/mainfig.gtex.raw.pdf'

# set parameters
nrows = 100; ncols = 100
highlightedtfs = [
        'FOXA1',
        'FOXA2',
        'EBF1',
        'HNF4G',
        'HNF4A',
        'PU1',
        'FOS',
        'TBP',
        'TAF1',
        'POL2',

        'GABPA',
        'CEBPB',
        'MAFF',
        'CTCF',
        'YY1'
        ]
highlightedtissues = [
        'Pancreas',
        'GI (liver)',
        'GI (stomach)',
        'Lymphocytes',
        'Spleen',
        'Fibroblasts',
        'Tibial nerve',

        'SM Artery (aorta)',
        'Brain (sub. nig.)',
        'Skeletal muscle',
        'Brain (putamen)',
        'SM Artery (tibial)'
        ]

# set up figure
fig = plt.figure(figsize=(4.06,7.1))
gs = gridspec.GridSpec(nrows,ncols)

# get data
## read results file and compute extreme p-values approximately if necessary
r = pd.read_csv(indir+'fdr5.gwresults', sep='\t')
r['id'] = [a.split(',')[0] for a in r.annot]
r['z_approx'] = r.mu / r['se(mu)']
lowp = (r.p<=1e-5)
r.loc[lowp, 'p'] = st.chi2.sf(r[lowp].z_approx**2, 1)
r['polp'] = -np.log10(r.p) * np.sign(r.z)
## read in tissue-specific results and merge in
resid = pd.read_csv(indir+'fdr5_tissuespecific.gwresults', sep='\t')
resid.pheno = resid.pheno.str.split('_resid').str.get(0)
resid['tissue_specific'] = 1
r = pd.merge(r, resid[['annot','pheno','tissue_specific']], how='left',
        on=['annot','pheno']).fillna(0)
## read in metadata and merge in
id2exp = pd.read_csv(params.sldp+'0.annotsummary/annotsummary.tsv', sep='\t')
exp2gene = pd.read_csv(params.sldp+'0.annotsummary/experiment_to_gene.tsv', sep='\t')
r = pd.merge(r, id2exp, on='id', how='left')
r = pd.merge(r, exp2gene, on='experiment', how='left')
r['phenoname'] = [info.phenotypes[p] for p in r.pheno]
## collapse different annots with same TF
toplot = r.groupby(['gene','phenoname']).agg(
        {
            'polp': lambda x : np.max(x) if np.max(x) > -np.min(x) else np.min(x),
            'tissue_specific': np.max
        }).reset_index()
print(toplot.shape)
print(toplot.tissue_specific.sum())
## create table for heatmap, then sort rows and columns by number of hits
tab = toplot.pivot(index='gene', columns='phenoname', values='polp').fillna(0)
tab['count'] = np.linalg.norm(tab.values, ord=0, axis=1)
tab['genename'] = tab.index
tab = tab.sort_values(by=['count', 'genename']).drop(['count', 'genename'], axis=1)
tab = tab.T
tab['count'] = np.linalg.norm(tab.values, ord=0, axis=1)
tab['phenoname'] = tab.index
tab = tab.sort_values(by=['count', 'phenoname']).drop(['count', 'phenoname'], axis=1)
## create binary table for tissue specific results
tab_resid = toplot.pivot(
        index='gene', columns='phenoname', values='tissue_specific').fillna(0).T
tab_resid = tab_resid.loc[tab.index][tab.columns]
## rotate tables to flip figure
tab = tab.T
tab_resid = tab_resid.T

## plot heatmap
from mpl_toolkits.axes_grid1 import make_axes_locatable
ax = plt.subplot(gs[:,:95])
# cax = plt.subplot(gs[95:,30:70])
cb = ax.matshow(tab.values, interpolation='nearest', vmin=-7, vmax=7,
        cmap='bwr', aspect='equal')
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cax.set_aspect(aspect=10)
## Major ticks
ax.set_xticks(range(len(tab.columns)))
ax.set_yticks(range(len(tab)))
## Labels for major ticks
ax.set_xticklabels(tab.columns, rotation=90, fontsize=5.5)
ax.set_yticklabels(tab.index, fontsize=5.5)
## Minor ticks
ax.set_xticks(np.arange(-.5, len(tab.columns), 1), minor=True);
ax.set_yticks(np.arange(-.5, len(tab), 1), minor=True);
## turn off the tick lines along each axis
ax.tick_params(
    axis='both',       # changes apply to both axes
    which='both',      # both major and minor ticks are affected
    pad=-4,
    left='off',
    right='off',
    bottom='off',      # ticks along the bottom edge are off
    top='off')         # ticks along the top edge are off
for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(0.1)
## Gridlines based on minor ticks
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.1)
## place asterisks on tissue-specific results
asterisks_y, asterisks_x = np.where(tab_resid)
ax.scatter(asterisks_x, asterisks_y, s=0.8, facecolors='none', edgecolors='black',
        linewidth=1.0)
## make highlighted relationships have bold labels
for l in ax.get_xticklabels() + ax.get_yticklabels():
    if l.get_text() in highlightedtfs+highlightedtissues:
        l.set_weight('extra bold')
## axis limits
ax.set_xlim(0-0.5, len(tab.columns)-0.5)
ax.set_ylim(len(tab.index)-0.5, 0-0.5)
## color bar
cax.tick_params(size=0, labelsize=7, pad=3)
cb = fig.colorbar(cb, cax=cax, ticks=[-7,0,7])
cb.ax.set_yticklabels([r'$-7$',r'$0$',r'$7$'])
cb.set_label(label=r'$-\log_{10}(p)$ (polarized)', fontsize=params.labelfontsize, labelpad=2)
cb.outline.set_linewidth(0.05)

# finish up (without despine)
fig.tight_layout()
print('saving figure', outname)
fs.makedir_for_file(outname)
plt.savefig(outname, dpi=500)
plt.close()
