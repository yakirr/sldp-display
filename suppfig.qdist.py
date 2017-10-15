from __future__ import print_function, division
import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import pyutils.fs as fs
import statutils.vis as vis
from plot import params, results_overview, results_detail, info

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/7.p9_a9/verboseresults/'
outname = me+'/out/suppfig.qdist.pdf'

# set parameters
toplot = results_overview.init(
        [indir+'../compiled_results/p9.complex.fdr5.indep'],
        [indir+'../compiled_results/p9.complex.fdr5.indep'],
        'gene').sort_values('sf_p')
toplot = toplot[['annot','pheno','gene', 'cell_line']].reset_index(drop=True)
nrows = 3; ncols = 4

# set up figure
fig = plt.figure(figsize=(6.5, 4.875))
gs = gridspec.GridSpec(nrows,ncols)

# make figure
for i,r in toplot.iterrows():
    # make plot of q
    ax = plt.subplot(gs[i])
    results_detail.plot_q(ax, indir, r.pheno, r.annot)
    ax.set_title('{} {} ({})'.format(info.phenotypes[r.pheno],
        r.gene, r.cell_line),
        fontsize=params.labelfontsize-1)

# finish up
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
