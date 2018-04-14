from __future__ import print_function, division
import os
import pandas as pd
from plot import params, info, results_overview, tables

me = os.path.dirname(os.path.abspath(__file__))
infile = params.sldprev + '1.basset1tfs_p12/old_results/complex'
outfile = me+'/out/supptable.complextraits_distincthits.tex'

# read in data
global results, passed
results = results_overview.init(
        [infile+'/fdr5_distinct.sub_ng.gwresults'],
        [infile+'/fdr5_distinct.sub_ng.gwresults'],
        'cell_line').sort_values('p')

# format
tables.format_sldp_results(results)
results['TF (num)'] = ['{} ({:d})'.format(g, n)
        for g,n in zip(results.gene, results['count'])]
results['Sign'] = ['+' if r>0 else r'\---' for r in results.rf]

results.enrichment = ['{:.0%}'.format(x).replace('%',r'\%') for x in results.enrichment]
results.cell_line = [x + ' ()' for x in results.cell_line]
results['Top 2 enrichments'] = ''

results = results[['Trait','TF (num)','rf','p','q']]

# print as tex
results.rename(
    columns={
        'p':'$p$',
        'rf':'$r_f$',
        'q':'$q$',
        }).to_latex(outfile,
                index=False,
                column_format='llrcc',
                bold_rows=True,
                escape=False)
