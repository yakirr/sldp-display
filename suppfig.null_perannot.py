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

sldp = '/groups/price/yakir/sldp/'
me = os.path.dirname(os.path.abspath(__file__))
indir = sldp+'/1.null_calib_a9/compiled_results'
outname = me+'/out/suppfig.null_perannot.png'

weights='Winv_ahat_h'

tickprops = {
        'direction':'out',
        'length':2,
        'width':0.8,
        'pad':4,
        'labelsize':6}
qqprops = {
        's':1.5,
        'fontsize':7,
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

## create qqplot for part a
# get data and write to file
print('creating part a')
sim = 'no_enrichment_3pcsaffy'
fname = '{}/{}.{}_{}.all'.format(
    indir, sim, 'maf5', weights)
print(fname)
processed = process(pd.read_csv(fname, sep='\t'))
processed.to_csv('{}/nullperannot.a.data.tsv'.format(indir), sep='\t')
print(sig.simes(processed.simes))
# make figure
vis.qqplot(processed.simes, errorbars=False, ax=ax1, **qqprops)
# make tick labels at integer increments
ax1.set_xticks(list(set(range(-5,10)) & set(ax1.get_xticks().astype(int))))
ax1.set_yticks(list(set(range(-5,10)) & set(ax1.get_yticks().astype(int))))
ax1.tick_params(**tickprops)

# create qqplot for part b
# get nobase data and write to file
print('creating part b')
sim = 'minor1_signed10_3pcsaffy'
fname = '{}/{}.{}_{}.all'.format(
    indir, sim, 'nobase', weights)
print(fname)
processed = process(pd.read_csv(fname, sep='\t'))
processed.to_csv('{}/nullperannot.b.nobase.data.tsv'.format(indir), sep='\t')
print(sig.simes(processed.simes))
# make figure
vis.qqplot(processed.simes, errorbars=False, c='r', label='no covariates', ax=ax2, **qqprops)
# get maf5 data and write to file
print('creating part b')
sim = 'minor1_signed10_3pcsaffy'
fname = '{}/{}.{}_{}.all'.format(
    indir, sim, 'maf5', weights)
print(fname)
processed = process(pd.read_csv(fname, sep='\t'))
processed.to_csv('{}/nullperannot.b.maf5.data.tsv'.format(indir), sep='\t')
print(sig.simes(processed.simes))
# make figure
vis.qqplot(processed.simes, errorbars=False, c='b', label='5 MAF bins', ax=ax2, **qqprops)
# make tick labels at integer increments
ax2.set_xticks(list(set(range(-5,10)) & set(ax2.get_xticks().astype(int))))
ax2.set_yticks(list(set(range(-5,10)) & set(ax2.get_yticks().astype(int))))
ax2.tick_params(**tickprops)
# legend
ax2.legend(loc='upper left', fontsize=5, markerscale=2, borderpad=0.1,
        labelspacing=0.2, columnspacing=0.2)
ax2.tick_params(**tickprops)

# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
# plt.show()
plt.savefig(outname, dpi=300); plt.close()
