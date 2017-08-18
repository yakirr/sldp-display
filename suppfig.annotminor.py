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
indir = params.sldp+'/0.annotsummary/'
outname = me+'/out/suppfig.annotminor.pdf'

# set aesthetics

# set params

## set up figure
fig = plt.figure(figsize=(3,3))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])

## get data
fname = '{}/{}'.format(
    indir, 'annotminor.tsv')
print(fname)
data = pd.read_csv(fname, sep='\t')

## create scatter plot
# get data
x = data['mean']; y = -np.log10(data.p)
err = data['stdev']
# scatter
ax1.errorbar(x, y, xerr=err, errorevery=1, fmt='o', mfc='b', ms=1.5, mew=0,
        linestyle='None',
        linewidth=0.5)

# vertical line at zero and FDR line (which is so low it doesnt show up)
passed = sig.sigrows(data, data.p, st=False, threshold=0.05)
print(len(passed))
ythresh = -np.log10(passed.p.max())
print(ythresh)
ax1.axhline(y=ythresh, **params.sig_thresh_line_props)
ax1.axvline(x=0, color='black', linewidth=0.5, alpha=0.8)

# axis limits and labels
ax1.set_xlim(-0.15, 0.15)
ax1.set_xticks([-0.15, 0, 0.15])
ax1.set_xlabel(r'$\bar{v}$ among non-zero entries', fontsize=params.labelfontsize)
ax1.set_ylabel(r'$-\log_{10}(p)$', fontsize=params.labelfontsize)
ax1.tick_params(**params.tickprops)


# finishing touches and save
sns.despine()
plt.tight_layout()

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname); plt.close()
