from __future__ import print_function, division
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from plot import params

def bias_plot(ax, indir, desc, weights, refpanel, estimand):
    # get data
    results = pd.read_csv(
            '{}{}.{}_{}.results'.format(indir, refpanel, desc, weights),
            sep='\t')
    print(refpanel, desc, weights, ':', len(results), 'results')

    # read in truth values
    truth = pd.read_csv(indir+'truth', sep='\t').set_index('simname')
    sims = truth.index.values

    # throw out some x-axis value to make the axis spacing right
    if estimand in ['h2v_h2g', 'h2v']:
        sims = np.concatenate([sims[0:1],sims[4:5],sims[6:]])

    # populate arrays to be plotted
    x = np.array([truth.loc[s,estimand] for s in sims])
    ests = np.array([results.loc[results.simname==s, estimand].values for s in sims]).T
    y = ests.mean(axis=0)
    stderr_mean = ests.std(axis=0) / np.sqrt(ests.shape[0])

    # do the plotting
    ax.boxplot(ests,
            positions=x,
            medianprops=dict(linewidth=0),
            widths=0.5*(np.max(x) - np.min(x))/len(x),
            labels=['{:.2g}'.format(x_) for x_ in x],
            boxprops={'linewidth':0.5},
            flierprops={'linewidth':0.5},
            capprops={'linewidth':0.5},
            whiskerprops={'linewidth':0.5, 'dashes':[2, 2]})
    ax.set_xticklabels([l.get_text() for l in ax.get_xticklabels()[::2]])
    ax.set_xticks(ax.get_xticks()[::2])
    ax.errorbar(x, y, yerr=stderr_mean,
            color='red',
            fmt='o',
            capthick=0.5,
            capsize=1.5,
            linewidth=0.5,
            markersize=1)

    # add y=x line
    marginx = 0.1*(max(x)-min(x))
    marginy = 0.1*(np.max(ests)-np.min(ests))
    ax.plot([0, max(x)+marginx], [0, max(x)+marginx],
            '--', color='black', label='truth', linewidth=0.5, dashes=[2,2])

    # format plot
    ax.axis((-marginx/2, max(x)+marginx,
            np.min(ests)-marginy/2, np.max(ests)+marginy))

def power_plot(ax, indir, desc, weights, refpanel,
        truthcolname, truthformat, labelfontsize=6, **powererrorbarprops):
    # get data
    results = pd.read_csv(
            '{}{}.{}_{}.results'.format(indir, refpanel, desc, weights),
            sep='\t')
    print(refpanel, desc, weights, ':', len(results), 'results')

    # add in true value of the parameter of interest
    truth = pd.read_csv(indir+'truth', sep='\t').set_index('simname')
    for s in truth.index.values:
        results.loc[results.simname==s, 'truth'] = \
                truthformat.format(truth.loc[s, truthcolname])

    # make power curve
    power = results[['truth']].drop_duplicates().reset_index()
    power['power'] = 0
    power['power_se'] = 0
    for i, row in power.iterrows():
        mask = (results['truth'] == row['truth'])
        p = power.loc[i,'power'] = (results[mask].sf_p <= 0.05).sum() / mask.sum()
        power.loc[i,'power_se'] = np.sqrt(p*(1-p)/mask.sum())
    power.sort_values(by='truth', inplace=True)

    # generate plot
    ax.errorbar(power['truth'], power.power,
            yerr=power.power_se,
            **powererrorbarprops)
    ax.set_ylabel('Power ($\\alpha=0.05$)', fontsize=labelfontsize)
