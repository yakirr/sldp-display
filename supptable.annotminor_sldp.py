from __future__ import print_function, division
import os
import pandas as pd
from plot import params, results_overview

me = os.path.dirname(os.path.abspath(__file__))
infile = params.sldp + '5.mafonly_a9/compiled_results'
outfile = me+'/out/supptable.annotminor_sldp.tex'

# read in data
global results, passed
results = results_overview.init(
        [infile+'/nobase_Winv_ahat_h.all'],
        [infile+'/nobase_Winv_ahat_h.fdr'])
passed = results[results.passed].copy()

# format
passed.r_f = ['{:.3f}'.format(x) for x in passed.r_f]
passed = passed[['origin','cell_line','gene','r_f','sf_p']]

# print as tex
passed.rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'gene':'Gene',
        'r_f':'$r_f$',
        'sf_p':'$p$'}).to_latex(outfile,
                index=False,
                bold_rows=True,
                escape=False)
