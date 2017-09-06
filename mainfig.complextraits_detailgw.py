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
from plot import params, results_detail

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/7.p9_a9/verboseresults/'
outname = me+'/out/mainfig.complextraits_detailgw.raw.pdf'

# set parameters
toplot = [
        # ('CD', 'SydhGm18951Pol2Iggmus', 'POL2'),
        # ('UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED','UtaMcf7CtcfSerumstim', 'CTCF'),
        ('UKB_460K.cov_EDU_YEARS','HaibGm12878Bcl11aPcr1x', 'BCL11A', 'LCL',
            (-1, 1, 1), (-4, 4, 4),
            (20, 21),
            2, 60.5, 61,
            60.684, 60.780,
            None), #rs10189857, at chr2:60713235, has chisq = 65.5015
        ('PASS_Lupus','SydhK562CtcfbIggrab', 'CTCF', 'K562',
            (-5, 5, 1), (-2, 2, 3),
            (20, 21),
            16, 67.096, 68.173,
            67.596, 67.673,
            None),
        ('CD', 'SydhK562Irf1Ifng6h', 'IRF1', 'K562',
            (-1,1,1), (-2, 2, 3),
            (20, 42),
            5, 131.3, 132.323,
            131.817, 131.826,
            '/groups/price/yakir/dpos/WEIGHTS/NTR.BLOOD.RNAARR/NTR.IRF1.wgt.RDat')
        ]
nrows = 3; ncols = 100
ahat_Rv_plot_ncols = 40

# set up figure
height_in=4
fig = plt.figure(figsize=(ncols/(ahat_Rv_plot_ncols*nrows)*height_in,height_in))
gs = gridspec.GridSpec(nrows,ncols)

# make figure
for i,(pheno, annot, tf, cell_line,
        (xmin, xmax, xexp), (ymin, ymax, yexp),
        (ytick, yrange),
        c, start, end, gstart, gend, twas) in enumerate(toplot):
    # make plot of ahat vs Rv
    results_detail.plot_ahat_vs_Rv(plt.subplot(gs[i,:ahat_Rv_plot_ncols]),
            indir, pheno, annot,
            (xmin, xmax, xexp),
            (ymin, ymax, yexp))

    # make manhattan plot
    results_detail.manhattan(fig, gs[i,ahat_Rv_plot_ncols:],
            pheno, tf,
            c, start, end,
            gstart, gend,
            ytick=ytick, yrange=yrange,
            twas=twas, show_gene_loc=False)

# finish up
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
