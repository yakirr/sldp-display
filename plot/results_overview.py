from __future__ import print_function, division
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import scipy.cluster.hierarchy as sch
import pandas as pd
import numpy as np
import info, params

volcanoprops = {
        'linewidth':0}
nonsig_color = (0.5, 0.5, 0.5, 0.4)

# convert long annot names to short annot names
def format_annot_names(df):
    df.annot = df.annot.str.split(',').str.get(0)
    return df

# add color metadata to annotations
def get_color(annots, colorby, apply_to=None):
    if apply_to is None:
        apply_to = [True]*len(annots)
    apply_to = np.array(apply_to, dtype=bool)

    colors = sns.hls_palette(
            n_colors=len(annots.loc[apply_to, colorby].unique()),
            )
    colordict = {a:nonsig_color for a in annots.loc[~apply_to, colorby].unique()}
    colordict.update(dict(zip(
        annots.loc[apply_to, colorby].unique(),
        [c+(0.8,) for c in colors])))
    return [colordict[x] for x in annots[colorby]]

# load up information for plot
def init(all_results, fdr_results, colorby):
    # read in results and merge in fdr information
    fdr = format_annot_names(
            pd.concat([pd.read_csv(f, sep='\t') for f in fdr_results], axis=0))
    fdr['passed'] = True
    fdrcols = list(set(['pheno','annot','passed','q'])&set(fdr.columns)) # in case q is absent
    results = pd.merge(
            format_annot_names(
                pd.concat([pd.read_csv(f, sep='\t') for f in all_results], axis=0)
                ).drop('q',axis=1, errors='ignore'),
            fdr[fdrcols],
            on=['pheno','annot'], how='left').fillna(False)

    # decide which annots to color and add color information to annot df
    annots_to_color = fdr.annot.unique()
    info.annots['color'] = get_color(info.annots, colorby,
            apply_to=[a in annots_to_color for a in info.annots.annot])

    # merge everything
    results = pd.merge(results, info.annots, on='annot', how='inner')
    phenos = results.loc[results.passed].pheno.unique()
    results = pd.merge(results, pd.DataFrame(phenos, columns=['pheno']),
            on='pheno', how='inner')
    print(len(results.pheno.unique()), 'phenos',
            len(results.annot.unique()), 'annots',
            len(results), 'results')
    return results

# single volcano plot
def volcano(ax, results, pheno, fontsize):
    # prepare data
    myresults = results[results.pheno == pheno].copy()
    myresults['logp'] = -np.log10(myresults.sf_p)
    myresults.loc[myresults.passed, 'markersize'] = 15
    myresults.loc[~myresults.passed, 'markersize'] = 5

    # gray out any TF with no passing annotations
    myresults['graymask'] = ~myresults.passed
    for gg in myresults[myresults.passed].genegroup.unique():
        myresults.loc[myresults.genegroup == gg, 'graymask'] = False
    myresults.loc[myresults.graymask, 'color'] = \
            myresults.loc[myresults.graymask, 'color'].apply(lambda x: nonsig_color)

    # sort so that TFs with fewer experiments that passed will get plotted on top
    counts = pd.DataFrame(
            myresults[myresults.passed].genegroup.value_counts()).rename(
                    columns={'genegroup':'count'})
    myresults = pd.merge(myresults,
            counts,
            left_on='genegroup', right_index=True, how='left').sort_values(
            ['graymask','count'], ascending=[False,False])

    # figure out where to draw FDR threshold
    high = myresults[myresults.passed].logp.min()
    low = myresults[~myresults.passed].logp.max()
    thresh = (high + low)/2

    # plot and create FDR line
    ax.scatter(myresults.r_f*100, myresults.logp,
            color=myresults.color, s=myresults.markersize,
            **volcanoprops)
    ax.axhline(y=thresh, **params.sig_thresh_line_props)

    # set labels
    ax.set_xlabel(r'Estimated $r_f$ (%)', fontsize=fontsize)
    ax.set_ylabel(r'$-\log_{10}(p)$', fontsize=fontsize)
    ax.set_title(info.phenotypes[pheno], fontsize=fontsize+1)

    # set axis limits and ticks
    xmin = 1.2*min(myresults.r_f*100); xmax = 1.2*max(myresults.r_f*100)
    delta = (xmax - xmin)/4
    ax.set_xticks(np.concatenate([
        np.arange(0, 1.2*xmin, -delta)[::-1], np.arange(delta, 1.2*xmax, delta)
        ]))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
    ax.set_xlim(xmin, xmax)
    ax.set_yticks(list(set(range(-5,10)) & set(ax.get_yticks().astype(int))))
    ax.set_ylim(-0.5, ax.get_ylim()[1])

# NOTE: this function assumes that init has been called already
def legend_contents(groupby):
    patches = []
    sort_annot = info.annots.sort_values(groupby)
    groups = sort_annot.loc[sort_annot[groupby].notnull(), groupby].unique()
    for group in groups:
        c = info.annots[info.annots[groupby] == group].color.unique()[0]
        if c != nonsig_color:
            patches.append(mpatches.Patch(color=c, label=group))
    # patches.append(mpatches.Patch(color=nonsig_color, label='No sig. results'))

    return patches

def segmented_bar(ax, passed, phenos, extra_mask, extra_color, title, fontsize,
        unmarked_color='b'):
    # filter to correct set of results
    mask = [p in phenos for p in passed.pheno]
    myresults = passed[mask].copy()
    # add the right order onto the annotation descriptions and sort
    myresults.desc = pd.Categorical(myresults.desc,
            info.category_order, ordered=True)
    myresults.sort_values(['desc', 'gene'], ascending=[True, True], inplace=True)
    # print basic info
    print(pd.DataFrame(myresults.gene.value_counts()).sort_index())
    print(len(myresults), 'total results')

    # make bar chart
    for i,(_, row) in enumerate(myresults.iterrows()):
        ax.barh(0, 1, 1, color=row.color, left=i, linewidth=0.2, edgecolor='white')
        if row[extra_mask]:
            ax.barh(1.1, 1, 0.2, color=extra_color, left=i, linewidth=0)
        else:
            ax.barh(1.1, 1, 0.2, color=unmarked_color, left=i, linewidth=0)
    ax.set_ylim(0,1.3)
    ax.set_xlim(0, len(myresults))

    # make ticks really small and erase axes
    ax.tick_params(
        length=0.5,
        pad=0.5,
        labelsize=0.5)
    ax.axis('off')

    # set title
    ax.set_title(title, fontsize=fontsize)

# ######### things below this line are unused ##########
# def summary_table(results, pheno):
#     summary = pd.DataFrame()
#     myresults = results[(results.pheno == pheno) & results.passed & (results.sf_p <= 1)]
#     for g in myresults.gene.unique():
#         thisbe = myresults[myresults.gene == g]
#         summary = summary.append({
#             'pheno':pheno,
#             'gene':g,
#             'num':len(thisbe),
#             'pmin':thisbe.sf_p.min(),
#             'pmax':thisbe.sf_p.max(),
#             'r_f':thisbe.r_f.mean(),
#             'cells':','.join(thisbe.cell_line.values)},
#             ignore_index=True)
#     summary.sort_values(
#             'pmin', ascending=True, inplace=True)
#     return summary[['gene','pheno','pmin','pmax','num','r_f','cells']]

# def plot_with_corr(results, pheno):
#     myresults = results[
#             (results.pheno == pheno) & results.passed & (results.sf_z**2>9)].set_index('annot')
#     sig = myresults.index.unique()
#     corr = pd.read_csv(info.sldp+'/6.annotcorr_a9/results/all.corr', sep='\t', index_col=0)
#     corr.index = corr.index.str.split(',').str.get(0)
#     corr.columns = corr.columns.str.split(',').str.get(0)
#     corr = corr.loc[sig][sig]

#     Y = sch.linkage(corr.values, method='centroid') # if the above lines aren't working
#     Z = sch.dendrogram(Y, orientation='right', no_plot=True)
#     ind = corr.index.values[Z['leaves']]

#     corr = corr.loc[ind][ind]
#     myresults = myresults.loc[ind]

#     fig = plt.figure(figsize=(6,6))
#     gs = gridspec.GridSpec(4,4)
#     ax1 = plt.subplot(gs[0:3,0:3])
#     ax2 = plt.subplot(gs[3,0:3])

#     print(corr.shape)
#     sns.heatmap(corr, square=True, xticklabels=False, yticklabels=1,
#             ax=ax1, cbar=False)
#     # ax2.scatter(range(len(myresults)), -np.log10(myresults.sf_p))
#     ax2.plot(-np.log10(myresults.sf_p))
#     ax2.set_xlim(-1, len(myresults))
#     sns.despine()
#     gs.tight_layout(fig)
#     plt.show()

# def heatmap(results):
#     myresults = results[results.passed & (results.sf_p <= 1e-3)].copy()
#     myresults['logp'] = -np.log10(myresults.sf_p)
#     tab = myresults.pivot(index='annot', columns='pheno', values='logp').fillna(0)

#     sig = tab.index.unique()
#     corr = pd.read_csv(info.sldp+'/6.annotcorr_a9/results/all.corr', sep='\t', index_col=0)
#     corr.index = corr.index.str.split(',').str.get(0)
#     corr.columns = corr.columns.str.split(',').str.get(0)
#     corr = corr.loc[sig][sig]

#     Y = sch.linkage(corr.values, method='centroid') # if the above lines aren't working
#     Z = sch.dendrogram(Y, orientation='right', no_plot=True)
#     ind = corr.index.values[Z['leaves']]

#     corr = corr.loc[ind][ind]
#     tab = tab.loc[ind][['BP_mono_gene_nor_combat_peer_10',
#         'BP_neut_gene_nor_combat_peer_10',
#         'geneexp_total_NTR',
#         'BP_mono_K27AC_log2rpm_peer_10',
#         'BP_neut_K27AC_log2rpm_peer_10',
#         'BP_mono_K4ME1_log2rpm_peer_10',
#         'BP_neut_K4ME1_log2rpm_peer_10',
#         'BP_tcel_K4ME1_log2rpm_peer_10']]

#     print(tab.shape)
#     sns.heatmap(tab, yticklabels=1)
#     plt.tight_layout()
#     plt.show()
