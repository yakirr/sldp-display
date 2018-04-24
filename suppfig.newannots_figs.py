from __future__ import print_function, division
import os
import numpy as np
import pandas as pd
import statutils.sig as sig
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import ypy.fs as fs
import seaborn as sns
from plot import params, results_overview; reload(results_overview)

me = os.path.dirname(os.path.abspath(__file__))

def zscore_plot(ax, title, oth_name, ba_allres, ba_fdrres, oth_allres, oth_fdrres):
    ## load basset and other sets of results
    ba = results_overview.init(
        ba_allres,
        ba_fdrres,
        'gene', force_all_phenos=True).sort_values('p').rename(columns={
            'passed':'sig_ba'})
    oth = results_overview.init(
        oth_allres,
        oth_fdrres,
        'gene', force_all_phenos=True).sort_values('p').rename(columns={
            'passed':'sig_oth'})
    ## merge and compute correlation
    cols = ['annot','pheno','z','p','mu']
    both = pd.merge(ba[cols+['sig_ba']], oth[cols+['sig_oth']],
            on=['annot','pheno'], how='inner')
    corr = np.corrcoef(both.z_x, both.z_y)[0,1]
    ## make figure
    mask = ~both.sig_ba
    ax.scatter(both[mask].z_x, both[mask].z_y, s=2, alpha=0.2, c='blue')
    ax.scatter(both[both.sig_ba].z_x, both[both.sig_ba].z_y, color='red', alpha=0.5, s=4)
    vmin = -5; vmax = 5
    # ax.axhline(y=1.96, xmin=vmin, xmax=vmax, **params.sig_thresh_line_props)
    # ax.axhline(y=-1.96, xmin=vmin, xmax=vmax, **params.sig_thresh_line_props)
    ax.plot([vmin, vmax], [vmin, vmax], **params.sig_thresh_line_props)
    ax.text(-4.25, 4, r'$r = {:.2g}$'.format(corr), fontsize=params.labelfontsize+1)
    ax.set_xlim(vmin, vmax); ax.set_ylim(vmin, vmax)
    ax.set_xlabel('Basset z-score', fontsize=params.labelfontsize)
    ax.set_ylabel(oth_name + ' z-score', fontsize=params.labelfontsize)
    ax.set_title(title, fontsize=params.labelfontsize+3)
    ax.tick_params(**params.tickprops)

def sldpcorr_plot(ax, ba_allres, oth_allres, corrs):
    ## load basset and other sets of results
    ba = pd.concat([pd.read_csv(n, sep='\t') for n in ba_allres])
    oth = pd.concat([pd.read_csv(n, sep='\t') for n in oth_allres])
    oth['annot'] = [a.replace('.motif','') for a in oth.annot] # for motif annotations
    ## merge
    cols = ['annot','pheno','z']
    both = pd.merge(ba[cols], oth[cols],
            on=['annot','pheno'], how='inner')
    ## compute sldpcorrs
    corrs = corrs.copy()
    for a in corrs.index:
        this = both[both.annot == a+'.R']
        corrs.loc[a, 'sldpcorr'] = np.corrcoef(this.z_x, this.z_y)[0,1]
    ## make figure
    ax.scatter(corrs.acorr, corrs.sldpcorr, s=8, alpha=0.5, linewidth=0, c='blue')
    vmin = min(corrs.acorr.min(), corrs.sldpcorr.min()) - 0.05
    vmax = 1
    ax.plot([-1,1], [-1,1], **params.sig_thresh_line_props)
    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.text(vmin+0.03, vmax-0.05,
            '{0:.1f}%'.format((corrs.sldpcorr >= corrs.acorr).sum() / len(corrs)*100),
            fontsize=params.labelfontsize+1)
    ax.text(vmax-0.2, vmin+0.03,
            '{0:.1f}%'.format((corrs.sldpcorr < corrs.acorr).sum() / len(corrs)*100),
            fontsize=params.labelfontsize+1)
    ax.set_xlabel('correlation across SNPs', fontsize=params.labelfontsize)
    ax.set_ylabel('correlation across phenotypes', fontsize=params.labelfontsize)
    ax.tick_params(**params.tickprops)


####################################
####################################
# deepsea
####################################
####################################
# auprcs
print('deepsea auprc')
## set up figure
outname = me+'/out/suppfig.deepsea_auprc.pdf'
fig = plt.figure(figsize=(3,3))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
## read in data
bmeta = pd.read_csv('../data/annot/basset/basset1tfs.metadata.tsv', sep='\t')
dmeta = pd.read_csv('../data/annot/deepsea/ds1tfs.metadata.tsv', sep='\t')
auprcs = pd.merge(bmeta, dmeta, on='id')
## make figure
ax1.scatter(auprcs.auprc, auprcs.AUPRC, s=4, alpha=0.5)
ax1.plot([0,1],[0,1], **params.sig_thresh_line_props)
ax1.set_xlim(0,1)
ax1.set_ylim(0,1)
ax1.axhline(y=0.3, **params.enrichment_thresh_line_props)
ax1.axvline(x=0.3, **params.enrichment_thresh_line_props)
ax1.set_xlabel('Basset AUPRC', fontsize=params.labelfontsize)
ax1.set_ylabel('Deepsea AUPRC', fontsize=params.labelfontsize)
ax1.tick_params(**params.tickprops)
## finish and save
sns.despine()
plt.tight_layout()
fs.makedir_for_file(outname)
plt.savefig(outname); plt.close()

# trait sets for z-scores and annotcorr
traitsets = [
        ('BLUEPRINT', 'DeepSEA',
            [params.sldprev + '/1.basset1tfs_p12/molecular_BP/all.gwresults'],
            [params.sldprev + '/1.basset1tfs_p12/molecular_BP/fdr5.gwresults'],
            [params.sldprev + '/3.deepsea_p12/basset1qc/molecular_BP/all.gwresults'],
            [params.sldprev + '/3.deepsea_p12/basset1qc/molecular_BP/fdr5.gwresults']
            ),
        ('NTR', 'DeepSEA',
            [params.sldprev + '/1.basset1tfs_p12/molecular_NTR/all.gwresults'],
            [params.sldprev + '/1.basset1tfs_p12/molecular_NTR/fdr5.gwresults'],
            [params.sldprev + '/3.deepsea_p12/basset1qc/molecular_NTR/all.gwresults'],
            [params.sldprev + '/3.deepsea_p12/basset1qc/molecular_NTR/fdr5.gwresults']
            ),
        ('GTEx', 'DeepSEA',
            [params.sldprev + '/1.basset1tfs_p12/old_results/molecular_gtexv7_tissues/all.gwresults'],
            [params.sldprev + '/1.basset1tfs_p12/old_results/molecular_gtexv7_tissues/fdr5.gwresults'],
            [params.sldprev + '/3.deepsea_p12/basset1qc/molecular_gtexv7_tissues/all.gwresults'],
            [params.sldprev + '/3.deepsea_p12/basset1qc/molecular_gtexv7_tissues/fdr5.gwresults']
            ),
        ('Complex traits', 'DeepSEA',
            [params.sldprev + '/1.basset1tfs_p12/old_results/complex/all.sub_ng.gwresults'],
            [params.sldprev + '/1.basset1tfs_p12/old_results/complex/fdr5.sub_ng.gwresults'],
            [params.sldprev + '/3.deepsea_p12/basset1qc/complex/all.gwresults'],
            [params.sldprev + '/3.deepsea_p12/basset1qc/complex/fdr5.gwresults']
            )
        ]

# comparison of z-scores
print('deepsea z-scores')
## set up figure
outname = me+'/out/suppfig.deepsea_zscores.png'; dpi=500
fig = plt.figure(figsize=(6,6))
gs = gridspec.GridSpec(2,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[1,0])
ax4 = plt.subplot(gs[1,1])
## make figure
for ax, record  in zip([ax1, ax2, ax3, ax4], traitsets):
    zscore_plot(ax, *record)
## finish and save
sns.despine()
plt.tight_layout()
fs.makedir_for_file(outname)
plt.savefig(outname, dpi=dpi); plt.close()

# annotcorr vs sldpcorr
print('deepsea sldpcorr')
## set up figure
outname = me+'/out/suppfig.deepsea_sldpcorr.pdf'
fig = plt.figure(figsize=(3,3))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
## make figure
corrs = pd.read_csv(params.sldprev+'/3.deepsea_p12/annotcorr.tsv', sep='\t',
        index_col='id')
sldpcorr_plot(ax1,
        sum([t[2] for t in traitsets], []),
        sum([t[4] for t in traitsets], []),
        corrs)
## finish and save
sns.despine()
plt.tight_layout()
fs.makedir_for_file(outname)
plt.savefig(outname); plt.close()

####################################
####################################
# hocomoco
####################################
####################################
# trait sets for z-scores and annotcorr
traitsets = [
        ('BLUEPRINT', 'HOCMOCO',
            [params.sldprev + '/1.basset1tfs_p12/molecular_BP/all.gwresults'],
            [params.sldprev + '/1.basset1tfs_p12/molecular_BP/fdr5.gwresults'],
            [params.sldprev + '/2.hocomotifwbasset_p12/expdiff_sum/molecular_BP/all.gwresults'],
            [params.sldprev + '/2.hocomotifwbasset_p12/expdiff_sum/molecular_BP/fdr5.gwresults']
            ),
        ('NTR', 'HOCMOCO',
            [params.sldprev + '/1.basset1tfs_p12/molecular_NTR/all.gwresults'],
            [params.sldprev + '/1.basset1tfs_p12/molecular_NTR/fdr5.gwresults'],
            [params.sldprev + '/2.hocomotifwbasset_p12/expdiff_sum/molecular_NTR/all.gwresults'],
            [params.sldprev + '/2.hocomotifwbasset_p12/expdiff_sum/molecular_NTR/fdr5.gwresults']
            ),
        ('GTEx', 'HOCMOCO',
            [params.sldprev + '/1.basset1tfs_p12/old_results/molecular_gtexv7_tissues/all.gwresults'],
            [params.sldprev + '/1.basset1tfs_p12/old_results/molecular_gtexv7_tissues/fdr5.gwresults'],
            [params.sldprev + '/2.hocomotifwbasset_p12/expdiff_sum/molecular_gtexv7_tissues/all.gwresults'],
            [params.sldprev + '/2.hocomotifwbasset_p12/expdiff_sum/molecular_gtexv7_tissues/fdr5.gwresults']
            ),
        ('Complex traits', 'HOCMOCO',
            [params.sldprev + '/1.basset1tfs_p12/old_results/complex/all.sub_ng.gwresults'],
            [params.sldprev + '/1.basset1tfs_p12/old_results/complex/fdr5.sub_ng.gwresults'],
            [params.sldprev + '/2.hocomotifwbasset_p12/expdiff_sum/complex/all.gwresults'],
            [params.sldprev + '/2.hocomotifwbasset_p12/expdiff_sum/complex/fdr5.gwresults']
            )
        ]
# comparison of z-scores
print('hocomoco z-scores')
## set up figure
outname = me+'/out/suppfig.motif_zscores.png'; dpi=500
fig = plt.figure(figsize=(6,6))
gs = gridspec.GridSpec(2,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[1,0])
ax4 = plt.subplot(gs[1,1])
## make figure
for ax, record  in zip([ax1, ax2, ax3, ax4], traitsets):
    zscore_plot(ax, *record)
## finish and save
sns.despine()
plt.tight_layout()
fs.makedir_for_file(outname)
plt.savefig(outname); plt.close()

# annotcorr vs sldpcorr
print('hocomoco sldpcorr')
## set up figure
outname = me+'/out/suppfig.motif_sldpcorr.pdf'
fig = plt.figure(figsize=(3,3))
gs = gridspec.GridSpec(1,1)
ax1 = plt.subplot(gs[0,0])
## make figure
corrs = pd.read_csv(params.sldprev+'/2.hocomotifwbasset_p12/annotcorr.tsv', sep='\t',
        index_col='id')
sldpcorr_plot(ax1,
        sum([t[2] for t in traitsets], []),
        sum([t[4] for t in traitsets], []),
        corrs)
## finish and save
sns.despine()
plt.tight_layout()
fs.makedir_for_file(outname)
plt.savefig(outname); plt.close()
