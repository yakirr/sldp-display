from __future__ import print_function, division
import sys, os
import numpy as np
import scipy.stats as st
import pandas as pd
import ypy.fs as fs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from plot import params, info; reload(info)

me = os.path.dirname(os.path.abspath(__file__))
indir_sldp = params.sldprev+'/1.basset1tfs_p12/old_results/molecular_gtexv7_tissues/'
gtextpm_dir = params.sldprev+'/1.basset1tfs_p12/tpm_gtex/'
gtextpm_file = gtextpm_dir + 'GTEx.v6tfs.tpm.tissuemedians.tsv'
gtextpmtissuenames_file = gtextpm_dir + 'tissue_to_pheno.tsv'
outname = me+'/out/suppfig.gtex_expression.pdf'


# set parameters
nrows = 1; ncols = 1

# set up figure
fig = plt.figure(figsize=(3,3))
gs = gridspec.GridSpec(nrows,ncols)
ax1 = plt.subplot(gs[0,0])

# get data
## read sldp results
r = pd.read_csv(indir_sldp+'all.gwresults', sep='\t')
r['annot'] = [a.split(',')[0] for a in r.annot]
passed = pd.read_csv(indir_sldp+'fdr5.gwresults', sep='\t')
passed['annot'] = [a.split(',')[0] for a in passed.annot]
passed['sig'] = True
r = pd.merge(r, passed[['annot','pheno','sig']], how='left',
        on=['annot','pheno']).fillna(False)
## read in metadata and merge in
r = pd.merge(r, info.annots, on='annot', how='left')
r = pd.merge(r, info.gene_to_ensembl, on='gene', how='left')
tissue_names = pd.read_csv(gtextpmtissuenames_file,
        sep='\t', header=None, names=['tissue', 'pheno'])
r = pd.merge(r.drop(['tissue'], axis=1), tissue_names, on='pheno', how='left')
## read in tpm data and merge in
tpm = pd.read_csv(gtextpm_file, sep='\t')
# r['phenoname'] = [info.phenotypes[p] for p in r.pheno]
r = pd.merge(r, tpm, on='ensgid', how='left')

# process data for plot
global summary
summary = pd.DataFrame()

for i, t in enumerate(r.tissue.unique()):
    print(t)
    df = r.loc[(r.tissue == t), ['annot','gene','sig',t]]
    print(len(df), df.sig.sum())
    df = df.loc[df[t].notnull()]
    print(len(df))
    df = df.drop_duplicates(subset=['gene','sig']) # count each gene only once
    df = df.groupby('gene').agg({
        'annot':(lambda x:x.iloc[0]),
        'sig':max,
        t:np.mean, # doesnt matter cuz theyre all identical
        })
    print(len(df), df.sig.sum())

    if df.sig.sum() > 0:
        thresh = 5

        summary = summary.append({
            'nsig':df.sig.sum(),
            'nnotsig':(~df.sig).sum(),
            'nsig_expr':(df.loc[df.sig, t] >= thresh).sum(),
            'nnotsig_expr':(df.loc[~df.sig, t] >= thresh).sum()},
            ignore_index=True)

# print summary statistics and p-value for expressed/not-expressed analysis
p1 = summary.nsig_expr.sum() / summary.nsig.sum()
p2 = summary.nnotsig_expr.sum() / summary.nnotsig.sum()
p = (summary.nsig_expr.sum()+summary.nnotsig_expr.sum())/ \
        (summary.nsig.sum() + summary.nnotsig.sum())
z = (p1-p2)/np.sqrt(p*(1-p)*(1/summary.nsig.sum()+1/summary.nnotsig.sum()))
p = st.norm.sf(z, 0, 1)
print(p1, p2, z, p)

# create plot for expressed/not expressed
global x, y
x = summary.nnotsig_expr / summary.nnotsig
y = summary.nsig_expr / summary.nsig
s = summary.nsig
rgba_colors = np.zeros((len(x),4))
rgba_colors[:,2] = 1.0
rgba_colors[:,3] = s / (2*max(s)) + 0.5

vmin=0.2; vmax=1.1
ax1.scatter(x, y, color=rgba_colors, s=7)
ax1.plot([vmin,vmax],[vmin,vmax], linestyle='--', c='gray')
ax1.set_xlim(vmin, vmax); ax1.set_ylim(vmin, vmax)
ax1.set_xlabel('percent non-significant TFs expressed', fontsize=params.labelfontsize)
ax1.set_ylabel('percent significant TFs expressed', fontsize=params.labelfontsize)
ax1.tick_params(**params.tickprops)
ax1.text(vmin+0.05, vmax-0.05,
        str(int((y>=x).sum() / len(x)*100))+'%', fontsize=params.labelfontsize+1)
ax1.text(vmax-0.1, vmin+0.05,
        str(int((x>y).sum() / len(x)*100))+'%', fontsize=params.labelfontsize+1)

# finishing touches and save
sns.despine()
gs.tight_layout(fig)

fs.makedir_for_file(outname)
plt.savefig(outname); plt.close()
plt.close()
