from __future__ import print_function, division
import os
import pandas as pd
from plot import params, info, results_overview, tables

me = os.path.dirname(os.path.abspath(__file__))
infile = params.sldp + '7.p9_a9/compiled_results'
outfile = me+'/out/maintable.complextraits_indephits.raw.tex'

# read in data
global results, passed
results = results_overview.init(
        [infile+'/p9.complex.fdr5.indep'],
        [infile+'/p9.complex.fdr5.indep'],
        'cell_line').sort_values('sf_p')

# format
tables.format_sldp_results(results)
for i, (ind, r) in enumerate(
        results[(results['count']>1)&(results.Trait != 'Eczema')].iterrows()):
    results.loc[ind, 'gene'] += ''.join(['*']*(i+1))
results['Top TF (num)'] = ['{} ({:d})'.format(g, n)
        for g,n in zip(results.gene, results['count'])]
results['Sign'] = ['+' if r>0 else r'\---' for r in results.r_f]

results.enrichment = ['{:.0%}'.format(x).replace('%',r'\%') for x in results.enrichment]
results.q = [tables.format_p(x) for x in results.q]
results.cell_line = [x + ' ()' for x in results.cell_line]

results = results[['Trait','Top TF (num)','cell_line','r_f','sf_p','q']]

# print as tex
results.rename(
    columns={
        'cell_line':'Top cell line',
        'sf_p':'$p$',
        'q':'$q$',
        'r_f':'$r_f$',
        }).to_latex(outfile,
                index=False,
                column_format='lllrccc',
                bold_rows=True,
                escape=False)
