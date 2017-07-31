from __future__ import print_function, division
import os
import pandas as pd
from plot import params

me = os.path.dirname(os.path.abspath(__file__))
infile = params.sldp + '0.annotsummary/annotsummary.tsv'
outfile = me+'/out/supptable.annotsummary.raw.tex'

# read in data
annotsummary = pd.read_csv(infile, sep='\t', index_col=0)

# format
annotsummary.supp = ['{:.0f}'.format(x) for x in annotsummary.supp]
annotsummary.supp_percent = ['{:.2f}'.format(x) for x in annotsummary.supp_percent]
annotsummary.norm2 = ['{:.2f}'.format(x) for x in annotsummary.norm2]

# print as tex
annotsummary.rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'experiment':'Experiment',
        'auprc':'BASSET AUPRC',
        'supp':'$|v|_0$',
        'supp_percent':'$|v|_0/M$ (\\%)',
        'norm2':'$|v|_2$'}).to_latex(outfile,
                index=False,
                bold_rows=True,
                escape=False,
                longtable=True)
