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
outname = me+'/out/mainfig.complextraits_newrels_detailgw.raw.pdf'
numerical_outname = me+'/out/xlsxtable.complextraits_newrels_detailgw.xlsx'
writer = pd.ExcelWriter(numerical_outname)

# set parameters
toplot = [
        ('UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED','UtaMcf7CtcfSerumstim', 'CTCF', 'MCF7',
            (-5, 5, 1), (-4, 4, 4),
            (20, 25),
            3, 186.939, 187.963,
            187.439, 187.463,
            None),
        ('PASS_Anorexia','HaibHepg2Sp1Pcr1x', 'SP1', 'HEPG2',
            (-5, 5, 1), (-4, 4, 4),
            (20, 25),
            12, 53.773979, 53.810226,
            53.273979, 54.310226,
            None),
        ('CD','SydhGm18951Pol2Iggmus', 'POL2', 'GM18951',
            (-5, 5, 1), (-4, 4, 4),
            (20, 25),
            19, 1.086578, 1.095391,
            0.586578, 1.595391,
            None)
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
    numbers = results_detail.plot_ahat_vs_Rv(plt.subplot(gs[i,:ahat_Rv_plot_ncols]),
            indir, pheno, annot,
            (xmin, xmax, xexp),
            (ymin, ymax, yexp))
    numbers.to_excel(
            writer, chr(65+i)+'.1 GWAS vs SLDP',
            index=False)

    # make manhattan plot
    numbers = results_detail.manhattan(fig, gs[i,ahat_Rv_plot_ncols:],
            pheno, tf,
            c, start, end,
            gstart, gend,
            ytick=ytick, yrange=yrange,
            twas=twas, show_gene_loc=False)
    numbers.to_excel(
            writer, chr(65+i)+'.2 Manhattan plot',
            index=False)

# finish up
sns.despine()
gs.tight_layout(fig)

print('saving figure and writing numerical results')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
writer.save()
