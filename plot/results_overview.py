from __future__ import print_function, division
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np
import info

results_overview_props = {
        'linewidth':0}

# convert long annot names to short annot names
def format_annot_names(df):
    df.annot = df.annot.str.split(',').str.get(0)
    return df

# add color metadata to annotations
def get_color(annots, colorby='bigexp', apply_to=None):
    if apply_to is None:
        apply_to = [True]*len(annots)
    apply_to = np.array(apply_to, dtype=bool)

    colors = sns.hls_palette(
            n_colors=len(annots.loc[apply_to, colorby].unique()),
            )
    colordict = {a:(0.5, 0.5, 0.5, 0.4) for a in annots.loc[~apply_to, colorby].unique()}
    colordict.update(dict(zip(
        annots.loc[apply_to, colorby].unique(),
        [c+(0.8,) for c in colors])))
    return [colordict[x] for x in annots[colorby]]

# load up information for plot
def init(all_results, fdr_results):
    # read in results and merge in fdr information
    fdr = format_annot_names(
            pd.concat([pd.read_csv(f, sep='\t') for f in fdr_results], axis=0))
    fdr['passed'] = True
    results = pd.merge(
            format_annot_names(
                pd.concat([pd.read_csv(f, sep='\t') for f in all_results], axis=0)),
            fdr[['pheno','annot','passed']],
            on=['pheno','annot'], how='left').fillna(False)

    # decide which annots to color and add color information to annot df
    annots_to_color = fdr.annot.unique()
    info.annots['color'] = get_color(info.annots,
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
    myresults['markersize'] = 5
    myresults.loc[myresults.passed, 'markersize'] = 15

    # sort so that TFs with fewer experiments that passed will get plotted on top
    counts = pd.DataFrame(
            myresults[myresults.passed].bigexp.value_counts()).rename(
                    columns={'bigexp':'count'})
    myresults = pd.merge(myresults,
            counts,
            left_on='bigexp', right_index=True, how='left').sort_values(
            'count', ascending=False)

    # plot and set text labels
    ax.scatter(myresults.r_f, myresults.logp,
            color=myresults.color, s=myresults.markersize,
            **results_overview_props)
    ax.set_xlabel(r'$\hat r_f$', fontsize=fontsize)
    ax.set_ylabel(r'$-\log_{10}(p)$', fontsize=fontsize)
    ax.set_title(info.phenotypes[pheno], fontsize=fontsize+1)

    # set axis limits and ticks
    xmin = 1.2*min(myresults.r_f); xmax = 1.2*max(myresults.r_f)
    ax.set_xticks(np.linspace(xmin, xmax, 4))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    ax.set_xlim(xmin, xmax)
    ax.set_yticks(list(set(range(-5,10)) & set(ax.get_yticks().astype(int))))
