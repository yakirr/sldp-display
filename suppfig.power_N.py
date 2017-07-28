from __future__ import print_function, division
import os
import numpy as np
import pandas as pd
import pyutils.fs as fs
import matplotlib.pyplot as plt
import seaborn as sns
from plot import params

sldp = '/groups/price/yakir/sldp/'
me = os.path.dirname(os.path.abspath(__file__))
indir = sldp+'/4.vary_N/compiled_results/'
outname = me+'/out/suppfig.power_N.pdf'

tickprops = {
        'direction':'out',
        'length':2,
        'width':0.8,
        'pad':4,
        'labelsize':6}

desc='maf5'
weights='Winv_ahat_h'
refpanel='KG3.wim9nm'

# set up figure
fig = plt.figure(figsize=(3,3))

# get data
results = pd.read_csv(
        '{}{}.{}_{}.results'.format(indir, refpanel, desc, weights),
        sep='\t')
print(refpanel, desc, weights, ':', len(results), 'results')

# add in true value of the parameter of interest
truth = pd.read_csv(indir+'truth', sep='\t').set_index('simname')
for s in truth.index.values:
    results.loc[results.simname==s, '$N$'] = '{:d}'.format(int(truth.loc[s,'N']))

# make power curve
power = results[['$N$']].drop_duplicates().reset_index()
power['power'] = 0
power['power_se'] = 0
for i, row in power.iterrows():
    mask = (results['$N$'] == row['$N$'])
    p = power.loc[i,'power'] = (results[mask].sf_p <= 0.05).sum() / mask.sum()
    power.loc[i,'power_se'] = np.sqrt(p*(1-p)/mask.sum())
power.sort_values(by='$N$', inplace=True)

# generate plot
plt.errorbar(power['$N$'], power.power,
        yerr=power.power_se,
        c='b', linewidth=1)
plt.axis((0, 50000, 0, 1))
plt.gca().set_xlabel('Sample size', fontsize=7)
plt.gca().set_ylabel('Power ($\\alpha=0.5$)', fontsize=7)
plt.title('Power as a function of $N$', fontsize=8)
plt.gca().tick_params(**tickprops)


# finishing touches and save
sns.despine()
plt.tight_layout()

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
