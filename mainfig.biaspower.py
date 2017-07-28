from __future__ import print_function, division
import numpy as np
import os
import pandas as pd
import pyutils.fs as fs
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from plot import params

sldp = '/groups/price/yakir/sldp/'
me = os.path.dirname(os.path.abspath(__file__))
indir = sldp+'/2.vary_h2v/compiled_results/'
outname = me+'/out/mainfig.biaspower.pdf'

tickprops = {
        'direction':'out',
        'length':2,
        'width':0.8,
        'pad':4,
        'labelsize':6}

desc='maf5'
refpanel='KG3.wim9nm'
estimands = {
        'mu':'$\mu$',
        'r_f':'$r_f$',
        'h2v':'$h^2_v$',
        'h2v_h2g':'$h^2_v/h^2_g$'
        }
weights = {
        'Winv_ahat_h':'non-trivial weights',
        'Winv_ahat_I':'trivial weights'
        }

## set up figure
fig = plt.figure(figsize=(5,2.5))
gs = gridspec.GridSpec(1,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])

## get truth information for x-axis labels
truth = pd.read_csv(indir+'truth', sep='\t').set_index('simname')
truth['r_f'] = truth.mu * np.sqrt(truth.v2) / np.sqrt(truth.h2g)
truth['h2v_h2g'] = truth.r_f**2
truth['h2v'] = truth.h2v_h2g * truth.h2g

## make bias plot for part a
print('making bias plot')
weight = 'Winv_ahat_h'
m = 'r_f'

# read in data
results = pd.read_csv(
        '{}{}.{}_{}.results'.format(indir, refpanel, desc, weight),
        sep='\t')
print(refpanel, desc, weight, ':', len(results), 'results')
sims = truth.index.values

# throw out some x-axis value to make the axis spacing right
if m in ['h2v_h2g', 'h2v']:
    sims = np.concatenate([sims[0:1],sims[4:5],sims[6:]])

# populate arrays to be plotted
x = np.array([truth.loc[s,m] for s in sims])
ests = np.array([results.loc[results.simname==s, m].values for s in sims]).T
y = ests.mean(axis=0)
stderr_mean = 1.00*ests.std(axis=0) / np.sqrt(ests.shape[0])

# do the plotting
ax1.boxplot(ests,
        positions=x,
        medianprops=dict(linewidth=0),
        widths=0.5*(np.max(x) - np.min(x))/len(x),
        labels=['{:.2g}'.format(x_) for x_ in x],
        boxprops={'linewidth':0.5},
        flierprops={'linewidth':0.5},
        capprops={'linewidth':0.5},
        whiskerprops={'linewidth':0.5, 'dashes':[2, 2]})
ax1.set_xticklabels([l.get_text() for l in ax1.get_xticklabels()[::2]])
ax1.set_xticks(ax1.get_xticks()[::2])
ax1.errorbar(x, y, yerr=stderr_mean,
        color='red',
        label='methodname',
        fmt='o',
        capthick=0.5,
        capsize=1.5,
        linewidth=0.5,
        markersize=1)

# add y=x line
marginx = 0.1*(max(x)-min(x))
marginy = 0.1*(np.max(ests)-np.min(ests))
ax1.plot([0, max(x)+marginx], [0, max(x)+marginx],
        '--', color='black', label='truth', linewidth=0.5, dashes=[2,2])

# format plot
ax1.axis((-marginx/2, max(x)+marginx,
        np.min(ests)-marginy/2, np.max(ests)+marginy))
ax1.tick_params(**tickprops)
ax1.set_xlabel(r'True '+estimands[m], fontsize=6)
ax1.set_ylabel(r'Estimated '+estimands[m], fontsize=6)


## create power plot for part b
print('making power plot')
desc='maf5'
m='r_f'
name, hue, (r1,w1), (r2,w2) = ('contrastweights', 'Weight scheme',
            ('KG3.wim9nm','Winv_ahat_h'),('KG3.wim9nm','Winv_ahat_I'))

# decide which values of the stratified variable correspond to red vs blue
if hue == 'Weight scheme':
    blue = weights[w1]; red = weights[w2]
elif hue == 'Ref panel':
    blue = r1; red = r2

# read in results files from the two conditions and merge them
results1 = pd.read_csv(
        '{}{}.{}_{}.results'.format(indir, r1, desc, w1),
        sep='\t')
results1['Weight scheme'] = weights[w1]
results1['Ref panel'] = r1
results2 = pd.read_csv(
        '{}{}.{}_{}.results'.format(indir, r2, desc, w2),
        sep='\t')
results2['Weight scheme'] = weights[w2]
results2['Ref panel'] = r2
results = pd.concat([results1, results2], axis=0).rename(
        columns={'sf_z':'z-score'})
print(r1, desc, w1, ':', len(results1), 'results')
print(r2, desc, w2, ':', len(results2), 'results')

# add in data about the true value of the parameter of interest
for s in truth.index.values:
    results.loc[results.simname==s, 'True $r_f$'] = '{:.2g}'.format(truth.loc[s,m])

# compute power estimates
power = results[['Weight scheme', 'Ref panel', 'True $r_f$'
    ]].drop_duplicates().reset_index()
power['power'] = 0
power['power_se'] = 0
for i, row in power.iterrows():
    mask = (results['Weight scheme'] == row['Weight scheme']) & \
            (results['Ref panel'] == row['Ref panel']) & \
            (results['True $r_f$'] == row['True $r_f$'])
    p = power.loc[i,'power'] = (results[mask].sf_p <= 0.05).sum() / mask.sum()
    power.loc[i,'power_se'] = np.sqrt(p*(1-p)/mask.sum())

# separate things back out into two series to be plotted
power1 = power[power[hue] == blue]
power2 = power[power[hue] == red]

# do the plotting
errorbarprops = {
        'linewidth':0.5,
        'markersize':1}
ax2.errorbar(power1['True $r_f$'], power1.power,
        yerr=1.00*power1.power_se,
        c='b',
        label=blue,
        **errorbarprops)
ax2.errorbar(power2['True $r_f$'], power2.power,
        yerr=1.00*power2.power_se,
        c='r',
        label=red,
        **errorbarprops)

# format plot
ax2.axis((-0.005, 0.055, 0, 1.05))
ax2.legend(loc='upper left', fontsize=5, markerscale=2, borderpad=0.1,
        labelspacing=0.2, columnspacing=0.2)
# ax2.set_xticks(ax2.get_xticks()[::2])
ax2.tick_params(**tickprops)
ax2.set_xlabel(r'True $r_f$', fontsize=6)
ax2.set_ylabel('Power ($\\alpha=0.05$)', fontsize=6)


# finishing touches and save
sns.despine()
gs.tight_layout(fig)

print(m, ': saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
