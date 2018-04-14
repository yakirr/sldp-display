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
from plot import params, results_overview, results_detail, info; reload(results_detail)

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldprev+'/1.basset1tfs_p12/old_results/complex/'
indirldblocks = params.sldprev+'/6.complex_msigdb/chunks/' # 1700 blocks
indirchunks = params.sldp+'/7.p9_a9/verboseresults/' # 300 blocks
outname = me+'/out/suppfig.qdist.pdf'
outname_numloci = me+'/out/xlsxtable.complextraits_numloci.xlsx'
writer_numloci = pd.ExcelWriter(outname_numloci)

# set parameters
toplot = results_overview.init(
        [indir+'/fdr5_indep.sub_ng.gwresults'],
        [indir+'/fdr5_indep.sub_ng.gwresults'],
        'gene').sort_values('p')
toplot = toplot[['annot','pheno','gene', 'cell_line', 'origin', 'rf']].reset_index(drop=True)
nrows = 3; ncols = 4

# set up figure
fig = plt.figure(figsize=(6.5, 4.875))
gs = gridspec.GridSpec(nrows,ncols)

# make figure
global numloci
numloci = pd.DataFrame()
for i,r in toplot.iterrows():
    # make plot of q
    ax = plt.subplot(gs[i])
    percentagree = results_detail.plot_q(ax, indirchunks, r.pheno, r.annot, np.sign(r.rf))
    numloci_ = results_detail.num_loci(indirldblocks, r.pheno, r.annot, np.sign(r.rf))
    print(numloci_, 'loci')
    print(percentagree, 'of ldblocks agree with genomewide trend')
    numloci = numloci.append({
        'Trait':info.phenotypes[r.pheno],
        'TF':r.gene,
        'Cell line':r.cell_line,
        'Lab':r.origin,
        '$r_f$':r.rf,
        'percent':percentagree,
        'Min. num. loci':numloci_}, ignore_index=True)
    ax.set_title('{} {} ({})'.format(info.phenotypes[r.pheno],
        r.gene, r.cell_line),
        fontsize=params.labelfontsize-1)

print(numloci.percent.mean(), 'agreement on avg')
numloci[['Trait','TF','Cell line','Lab','$r_f$','Min. num. loci']
        ].to_excel(writer_numloci, 'Min. num. loci',
                index=False)
# finish up
sns.despine()
gs.tight_layout(fig)

print('saving figure')
# writer_numloci.save() # this info is in Table 1
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
