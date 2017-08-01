from __future__ import print_function, division
import sys, os
import numpy as np
import pandas as pd
import pyutils.fs as fs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from statutils import vis, sig
import seaborn as sns
from plot import params

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/1.null_calib_a9/compiled_results'
outname = me+'/out/mainfig.null_aggregate.raw.png'

# set aesthetics
qqprops = {
        's':1.5,
        'fontsize':params.labelfontsize,
        'linewidth':0}

# set params
weights='Winv_ahat_h'

## set up figure
fig = plt.figure(figsize=(6,2))
gs = gridspec.GridSpec(1,3)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[0,2])

## create qqplot for part a
print('creating part a')
sim = 'no_enrichment_3pcsaffy'
fname = '{}/{}.{}_{}.all'.format(
    indir, sim, 'maf5', weights)
print(fname)
results = pd.read_csv(fname, sep='\t')
print((results.sf_p <= 0.05).sum()/float(len(results)), sig.chi2(results.sf_p))
vis.qqplot(results.sf_p, errorbars=False, ax=ax1, **qqprops)
ax1.text(2.8, -0.5, 'avg $\chi^2$: {:.3g}'.format((results.sf_z**2).mean()),
        fontsize=params.labelfontsize, color='b')
# make tick labels at integer increments
ax1.set_xticks(list(set(range(-5,10)) & set(ax1.get_xticks().astype(int))))
ax1.set_yticks(list(set(range(-5,10)) & set(ax1.get_yticks().astype(int))))
ax1.tick_params(**params.tickprops)

## create qqplot for part b
print('creating part b')
sim = 'unsigned'
fname = '{}/{}.{}_{}.all'.format(
    indir, sim, 'maf5', weights)
print(fname)
results = pd.read_csv(fname, sep='\t', nrows=427*3)
print((results.sf_p <= 0.05).sum()/float(len(results)), sig.chi2(results.sf_p))
vis.qqplot(results.sf_p, errorbars=False, ax=ax2, **qqprops)
ax2.text(1.5, -0.25, 'avg $\chi^2$: {:.3g}'.format((results.sf_z**2).mean()),
        fontsize=params.labelfontsize, color='b')
# make tick labels at integer increments
ax2.set_xticks(list(set(range(-5,10)) & set(ax2.get_xticks().astype(int))))
ax2.set_yticks(list(set(range(-5,10)) & set(ax2.get_yticks().astype(int))))
ax2.tick_params(**params.tickprops)

# # create qqplot for part c
print('creating part c')
sim = 'minor1_signed10_3pcsaffy'
fname = '{}/{}.{}_{}.all'.format(
    indir, sim, 'nobase', weights)
print(fname)
results = pd.read_csv(fname, sep='\t')
print((results.sf_p <= 0.05).sum()/float(len(results)), sig.chi2(results.sf_p))
vis.qqplot(results.sf_p, errorbars=False, ax=ax3, c='r', label='no covariates', **qqprops)
ax3.text(2.8, 0.2, 'avg $\chi^2$: {:.3g}'.format((results.sf_z**2).mean()),
        fontsize=params.labelfontsize, color='r')

fname = '{}/{}.{}_{}.all'.format(
    indir, sim, 'maf5', weights)
print(fname)
results = pd.read_csv(fname, sep='\t')
print((results.sf_p <= 0.05).sum()/float(len(results)), sig.chi2(results.sf_p))
vis.qqplot(results.sf_p, errorbars=False, ax=ax3, label='5 MAF bins',  **qqprops)
ax3.text(2.8, -0.5, 'avg $\chi^2$: {:.3g}'.format((results.sf_z**2).mean()),
        fontsize=params.labelfontsize, color='b')
# add legend and format ticks
ax3.legend(loc='upper left', fontsize=5, markerscale=2, borderpad=0.1,
        labelspacing=0.2, columnspacing=0.2)
ax3.tick_params(**params.tickprops)
# make tick labels at integer increments
ax2.set_xticks(list(set(range(-5,10)) & set(ax2.get_xticks().astype(int))))
ax2.set_yticks(list(set(range(-5,10)) & set(ax2.get_yticks().astype(int))))

# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
# plt.show()
plt.savefig(outname, dpi=300); plt.close()
