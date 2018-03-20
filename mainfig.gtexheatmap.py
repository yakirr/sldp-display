from __future__ import print_function, division
import pandas as pd
import seaborn as sns
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats as st
from plot import params, results_overview

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/7.p9_a9/compiled_results/'
outname = me+'/out/mainfig.complextraits_overview.raw.pdf'

resultsfile = 'nonresid.gwresults'
residfile = 'all.gwresults'
pheno2system = pd.read_csv('../meta/gtex_tissues.tsv', header=None,
        sep='\t', names=['pheno','system'])
id2exp = pd.read_csv('/n/groups/price/yakir/sldp.sub_ng/0.annotsummary/annotsummary.tsv',
        sep='\t')
exp2gene = pd.read_csv('/n/groups/price/yakir/sldp.sub_ng/0.annotsummary/experiment_to_gene.tsv',
        sep='\t')
gtexnames = pd.read_csv('../meta/gtex_names.tsv',
        sep='\t', header=None, names=['pheno', 'name'], index_col=0)

global r
r = pd.read_csv(resultsfile, sep='\t')
r['id'] = [a.split(',')[0] for a in r.annot]
r.z = r.mu / r['se(mu)']
r.loc[r.p<=1e-5, 'p'] = st.chi2.sf(r[r.p<=1e-5].z**2, 1)
r['polp'] = -np.log10(r.p) * np.sign(r.z)
r['sign'] = np.sign(r.z)
r = pd.merge(r, id2exp, on='id', how='left')
r = pd.merge(r, exp2gene, on='experiment', how='left')
r = pd.merge(r, pheno2system, on='pheno', how='left')

resid = pd.read_csv(residfile, sep='\t')
resid.pheno = resid.pheno.str.split('_resid').str.get(0)
resid.z = resid.mu / resid['se(mu)']
resid.loc[resid.p<=1e-5, 'p'] = st.chi2.sf(resid[resid.p<=1e-5].z**2, 1)
r = pd.merge(r,
        resid[['annot','pheno','z','p']].rename(columns={'z':'z_resid','p':'p_resid'}),
        how='left',
        on=['annot','pheno'])
r['polp_resid'] = -np.log10(st.chi2.sf(r.z_resid**2, 1)) * np.sign(r.z_resid)

global passed, toplot
passed = sig.sigrows_strat(r, r.p, r.pheno)
print('before thresh', len(passed))
passed = passed[passed.p <= 0.05]
print('after thresh', len(passed))
# passed['asterisk'] = 0
# passed.loc[passed.p_resid < 1e-2, 'asterisk'] = 1
# passed.loc[np.sign(passed.z) != np.sign(passed.z_resid), 'asterisk'] = 0
# passed_resid = sig.sigrows_strat(r, r.p_resid, r.pheno)
# passed_resid = r[r.p_resid <= 0.05/382].copy()
passed_resid = passed[-np.log10(passed.p) - -np.log10(passed.p_resid) < 0].copy()
passed_resid['asterisk'] = 1
passed = pd.merge(passed, passed_resid[['annot','pheno','asterisk']], on=['annot','pheno'],
        how='left').fillna(0)
passed['tophit'] = 0
for p in passed.pheno.unique():
    passed.loc[np.argmin(passed[passed.pheno == p].p), 'tophit'] = 1

toplot = passed.groupby(['gene','pheno']).agg(
        {
            'polp': lambda x : np.max(x) if np.max(x) > -np.min(x) else np.min(x),
            'polp_resid': lambda x : np.max(x) if np.max(x) > -np.min(x) else np.min(x),
            'asterisk': np.max,
            'tophit': np.max
        }).reset_index()
print(toplot.shape)
print(toplot.asterisk.sum())

global tab, tab_resid, tab_tophit
tab = toplot.pivot(index='gene', columns='pheno', values='polp').fillna(0)
tab['count'] = np.linalg.norm(tab.values, ord=0, axis=1)
tab['genename'] = tab.index
tab = tab.sort_values(by=['count', 'genename']).drop(['count', 'genename'], axis=1)
tab = tab.T
tab['count'] = np.linalg.norm(tab.values, ord=0, axis=1)
tab['phenoname'] = tab.index
tab = tab.sort_values(by=['count', 'phenoname']).drop(['count', 'phenoname'], axis=1)

tab_resid = toplot.pivot(index='gene', columns='pheno', values='asterisk').fillna(0).T
tab_resid = tab_resid.loc[tab.index][tab.columns]
tab_tophit = toplot.pivot(index='gene', columns='pheno', values='tophit').fillna(0).T
tab_tophit = tab_tophit.loc[tab.index][tab.columns]

tab.index = gtexnames.loc[tab.index, 'name']


# set up figure
fig = plt.figure(figsize=(5,3))
gs = gridspec.GridSpec(10,100)

# plot matrix
ax = plt.subplot(gs[:,:90])
cax = plt.subplot(gs[2:8,95:98])
cb = ax.matshow(tab.values, interpolation='nearest', vmin=-7, vmax=7, cmap='bwr', aspect='equal')
# turn off the actual tick lines
ax.tick_params(
    axis='both',          # changes apply to both axes
    which='both',      # both major and minor ticks are affected
    left='off',
    right='off',
    bottom='off',      # ticks along the bottom edge are off
    top='off')         # ticks along the top edge are off
for axis in ['top','bottom','left','right']:
      ax.spines[axis].set_linewidth(0.1)
# Major ticks
ax.set_xticks(range(len(tab.columns)))
ax.set_yticks(range(len(tab)))
# Labels for major ticks
ax.set_xticklabels(tab.columns, rotation=90, fontsize=3)
ax.set_yticklabels(tab.index, fontsize=3)
# Minor ticks
ax.set_xticks(np.arange(-.5, len(tab.columns), 1), minor=True);
ax.set_yticks(np.arange(-.5, len(tab), 1), minor=True);
# Gridlines based on minor ticks
ax.grid(which='minor', color='gray', linestyle='-', linewidth=0.1)

# place asterisks on tissue-specific results
# asterisks_y, asterisks_x = np.where(tab_resid)
# ax.scatter(asterisks_x, asterisks_y, s=0.2, facecolors='none', edgecolors='black', linewidth=0.5)

# place circles on top hits within each phenotype
# circles_y, circles_x = np.where(tab_tophit)
# ax.scatter(circles_x, circles_y, s=4, facecolors='none', edgecolors='black', linewidth=0.2)
ax.set_xlim(0-0.5, len(tab.columns)-0.5)
ax.set_ylim(len(tab.index)-0.5, 0-0.5)

cax.tick_params(size=0, labelsize=4)
cb = fig.colorbar(cb, cax=cax)
cb.set_label(label=r'$-\log_{10}p$ (polarized)', fontsize=5)
cb.outline.set_linewidth(0.05)

# gs.tight_layout(fig)

# plt.show()
plt.savefig('out/heatmap.png', dpi=300)
plt.close()
