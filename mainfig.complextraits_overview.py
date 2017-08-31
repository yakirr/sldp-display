from __future__ import print_function, division
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyutils.fs as fs
from plot import params, results_overview; reload(results_overview)

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/7.p9_a9/compiled_results/'
outname = me+'/out/mainfig.complextraits_overview.raw.pdf'

# set parameters
phenos = [
    'CD',
    'PASS_Lupus',
    'UKB_460K.cov_EDU_YEARS',
    'UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED',
    # 'PASS_HDL',
    'PASS_Anorexia'
    ]
nrows = 2; ncols = 3

# set up figure
fig = plt.figure(figsize=(6,4))
gs = gridspec.GridSpec(nrows,ncols)

# get data
results = results_overview.init(
        [indir+'/p9.complex.all'],
        [indir+'/p9.complex.fdr5'],
        'genegroup')

# make figure
for cell, pheno in zip(gs, phenos):
    ax = plt.subplot(cell)
    print(pheno)
    results_overview.volcano(ax, results, pheno, params.labelfontsize)
    ax.tick_params(**params.tickprops)

# add legend in bottom right
ax = plt.subplot(gs[1,2])
patches = results_overview.legend_contents('genegroup')
ax.set_axis_off()
ax.legend(handles=patches, fontsize=6, markerscale=2, borderpad=0.1,
    labelspacing=0.4, columnspacing=0.2, loc='upper left')

# finish up
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
