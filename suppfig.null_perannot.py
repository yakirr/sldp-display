from __future__ import print_function, division
import sys, os
import numpy as np
import scipy.stats as st
import pandas as pd
import pyutils.fs as fs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from statutils import vis, sig
import seaborn as sns
from plot import params

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/1.null_calib_a9/compiled_results'
outname = me+'/out/suppfig.null_perannot.raw.pdf'

weights='Winv_ahat_h'

fontsize=7
tickprops = {
        'direction':'out',
        'length':2,
        'width':0.8,
        'pad':4,
        'labelsize':6}
qqprops = {
        's':1.5,
        'fontsize':fontsize,
        'linewidth':0}

def process(results):
    processed = pd.DataFrame(columns=['annot']).set_index('annot')
    for a in results.annot.unique():
        p = results.loc[results.annot==a, 'sf_p']
        processed.loc[a, 'simes'] = sig.simes(p)
        processed.loc[a, 'avgchi2'] = st.chi2.isf(p, 1).mean()
    processed.simes = pd.to_numeric(processed.simes)
    return processed

## set up figure
fig = plt.figure(figsize=(4,2))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])

## get data
print('creating part a')
sim = 'no_enrichment_3pcsaffy'
fname = '{}/{}.{}_{}.all'.format(
    indir, sim, 'maf5', weights)
print(fname)
processed = process(pd.read_csv(fname, sep='\t'))
mafres = pd.read_csv(params.sldp+'/5.mafonly_a9/compiled_results/nobase_Winv_ahat_h.all',
        sep='\t').rename(columns={'sf_z':'maf_z'})
processed = pd.merge(processed, mafres[['annot','maf_z']], left_index=True, right_on='annot',
        how='left')
processed.to_csv('{}/nullperannot.data.tsv'.format(indir), sep='\t')
print(sig.simes(processed.simes))

## create qqplot for part a
vis.qqplot(processed.simes, errorbars=False, ax=ax1, **qqprops)
# make tick labels at integer increments
ax1.set_xticks(list(set(range(-5,10)) & set(ax1.get_xticks().astype(int))))
ax1.set_yticks(list(set(range(-5,10)) & set(ax1.get_yticks().astype(int))))
ax1.tick_params(**tickprops)

## create scatter plot for part b
# get data
x = np.abs(processed.maf_z); y = processed.avgchi2
# scatter
ax2.scatter(x, y, s=1.5, linewidth=0)
# add trendline
ax2.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),
        c='gray', linewidth=0.8, dashes=[2,2])
r2 = np.corrcoef(processed.maf_z**2, processed.avgchi2)[0,1]**2
ax2.text(3, 0.82, '$r^2$: {:.2g}'.format(r2),
        fontsize=6, color='b')

ax2.set_xlabel(r'$|z|$, minor-allele-only trait', fontsize=fontsize)
ax2.set_ylabel(r'avg $\chi^2$', fontsize=fontsize)
ax2.tick_params(**tickprops)
print(r2)


# finishing touches and save
sns.despine()
plt.tight_layout()

print('saving figure')
fs.makedir_for_file(outname)
# plt.show()
plt.savefig(outname); plt.close()
