from __future__ import print_function, division
import os
import pandas as pd
from plot import params, info, results_overview

me = os.path.dirname(os.path.abspath(__file__))
infile = params.sldp + '7.p9_a9/compiled_results'
outfile = me+'/out/maintable.complextraits_indephits.tex'

# read in data
global results, passed
results = results_overview.init(
        [infile+'/p9.complex.fdr5.indep'],
        [infile+'/p9.complex.fdr5.indep'],
        'cell_line').sort_values('sf_p')

# format
results['Trait'] = [info.phenotypes[p] for p in results.pheno]
for i, (ind, r) in enumerate(
        results[(results['count']>1)&(results.Trait != 'Eczema')].iterrows()):
    results.loc[ind, 'gene'] += ''.join(['*']*(i+1))
results['TF (num)'] = ['{} ({:d})'.format(g, n)
        for g,n in zip(results.gene, results['count'])]
results['Sign'] = ['+' if r>0 else r'\---' for r in results.r_f]

def format_p(p):
    temp = '{:.1e}'.format(p)
    base = temp[:temp.index('e')]
    exp = str(int(temp[temp.index('e')+1:]))
    return '$'+base+'\\times 10^{'+exp+'}$'
results.sf_p = [format_p(x) for x in results.sf_p]
results.enrichment = ['{:.0%}'.format(x).replace('%',r'\%') for x in results.enrichment]
results.r_f = ['{:.1%}'.format(x).replace('%',r'\%') for x in results.r_f]
results.q = [format_p(x) for x in results.q]
results.cell_line = [x + ' ()' for x in results.cell_line]
results = results[['Trait','TF (num)','cell_line','r_f','sf_p','q']]

# print as tex
results.rename(
    columns={
        'cell_line':'Top cell line',
        'sf_p':'$p$',
        'q':'$q$',
        'r_f':'$r_f$',
        }).to_latex(outfile,
                index=False,
                column_format='lllcccc',
                bold_rows=True,
                escape=False)
