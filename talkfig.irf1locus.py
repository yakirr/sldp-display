from __future__ import print_function, division
import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import ypy.fs as fs
import statutils.vis as vis
from plot import params, results_detail

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/7.p9_a9/verboseresults/'
outname = me+'/out/talkfig.irf1locus.pdf'
print(outname)

# set parameters
toplot = [
        ('CD', 'SydhK562Irf1Ifng6h', 'IRF1', 'K562',
            (-1,1,1), (-2, 2, 3),
            (20, 42),
            5, 131.3, 132.323,
            131.817, 131.826,
            '/n/groups/price/yakir/dpos/WEIGHTS/NTR.BLOOD.RNAARR/NTR.IRF1.wgt.RDat')
        ]
ahat_Rv_plot_ncols = 40

# set up figure
fig = plt.figure(figsize=(8,3))
gs = gridspec.GridSpec(1,1)

# make figure
for i,(pheno, annot, tf, cell_line,
        (xmin, xmax, xexp), (ymin, ymax, yexp),
        (ytick, yrange),
        c, start, end, gstart, gend, twas) in enumerate(toplot):
    # make manhattan plot
    numbers = results_detail.manhattan(fig, gs[i,:],
            pheno, tf,
            c, start, end,
            gstart, gend,
            ytick=ytick, yrange=yrange,
            twas=twas, show_gene_loc=False)

# finish up
sns.despine()
gs.tight_layout(fig)

print('saving figure and writing numerical results')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
