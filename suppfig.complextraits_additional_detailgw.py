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
from plot import params, results_detail; reload(results_detail)

# set up latex text handling. NOTE: if running on O2, remove this code
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/7.p9_a9/verboseresults/'
enrichments_file = params.sldprev+'/6.complex_msigdb/fdr5.enrichments'
genes_file = '../data/reference/genes/master_autosomes.tsv'
outname = me+'/out/suppfig.complextraits_additional_detailgw.raw.pdf'
numerical_outname = me+'/out/xlsxtable.complextraits_additional_detailgw.xlsx'
writer = pd.ExcelWriter(numerical_outname)

# set parameters
toplot = [
        ('PASS_Lupus','SydhK562CtcfbIggrab', 'CTCF', 'K562',
            (-5, 5, 1), (-2, 2, 3),
            (20, 21),
            16, 67.096, 68.173,
            67.596, 67.673,
            None,
            ['Ctcf','Rad21']),
        ('CD','SydhGm18951Pol2Iggmus', 'POL2', 'GM18951',
            (-1.2, 1.2, 1), (-2.2, 2.2, 3),
            (20, 21),
            19, 0.586578, 1.595391,
            1.086578, 1.095391,
            None,
            ['Pol2','Tbp','Taf1']),
        ]
nrows = 2; ncols = 100
ahat_Rv_plot_ncols = 26
manhattan_ncols = int(1.5*ahat_Rv_plot_ncols)
enrichment_ncols = 30

# set up figure
height_in=(2/3)*4
width_in=ncols/(ahat_Rv_plot_ncols*nrows)*height_in
print('width =', width_in)
fig = plt.figure(figsize=(width_in,height_in))
gs = gridspec.GridSpec(nrows,ncols)

# make figure
genes = pd.read_csv(genes_file, sep='\t')
for i,(pheno, annot, tf, cell_line,
        (xmin, xmax, xexp), (ymin, ymax, yexp),
        (ytick, yrange),
        c, start, end, gstart, gend, twas,
        enrichment_tfs) in enumerate(toplot):
    # figure out number of genes in locus
    genes_in_locus = genes[
            (genes['Chromosome/scaffold name'] == c) &
            (genes['Gene start (bp)'] <= end*1e6) &
            (genes['Gene end (bp)'] >= start*1e6)]
    print(genes_in_locus[['Gene name','Gene type','Gene start (bp)', 'Gene end (bp)']])
    print(len(genes_in_locus[genes_in_locus['Gene type'] == 'protein_coding']),
        'protein coding genes in locus')

    # make plot of ahat vs Rv
    numbers = results_detail.plot_ahat_vs_Rv(plt.subplot(gs[i,:ahat_Rv_plot_ncols]),
            indir, pheno, annot,
            (xmin, xmax, xexp),
            (ymin, ymax, yexp))
    numbers.to_excel(
            writer, chr(65+i)+'.1 GWAS vs SLDP',
            index=False)

    # make manhattan plot
    numbers = results_detail.manhattan(fig,
            gs[i,ahat_Rv_plot_ncols:(ahat_Rv_plot_ncols+manhattan_ncols)],
            pheno, tf,
            c, start, end,
            gstart, gend,
            ytick=ytick, yrange=yrange,
            twas=twas, show_gene_loc=False)
    numbers.to_excel(
            writer, chr(65+i)+'.2 Manhattan plot',
            index=False)

    # make enrichment plot
    numbers = results_detail.enrichment(
            plt.subplot(gs[i, -enrichment_ncols:]),
            enrichments_file,
            pheno, enrichment_tfs,
            [1, 15, 30], 30)
    numbers.to_excel(
            writer, chr(65+i)+'.3 Top enrichments',
            index=False)

# finish up
sns.despine()
gs.tight_layout(fig)

print('saving figure and writing numerical results')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
writer.save()
