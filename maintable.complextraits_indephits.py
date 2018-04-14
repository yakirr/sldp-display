from __future__ import print_function, division
import os
import pandas as pd
import numpy as np
from plot import params, info, results_overview, tables, results_detail

me = os.path.dirname(os.path.abspath(__file__))
infile = params.sldprev + '1.basset1tfs_p12/old_results/complex'
indirldblocks = params.sldprev+'/6.complex_msigdb/chunks/' # 1700 blocks
outfile = me+'/out/maintable.complextraits_indephits.raw.tex'
outfile_xlsx = me+'/out/maintable.complextraits_indephits.raw.xlsx'
writer_xlsx = pd.ExcelWriter(outfile_xlsx)

# read in data and add number of loci to each result
global results, passed
results = results_overview.init(
        [infile+'/fdr5_indep.sub_ng.gwresults'],
        [infile+'/fdr5_indep.sub_ng.gwresults'],
        'cell_line').sort_values('p')
results[r'\# sites'] = [
        int(results_detail.num_loci(indirldblocks, r.pheno, r.annot, np.sign(r.rf)))
        for _, r in results.iterrows()]

# format
tables.format_sldp_results(results)
for i, (ind, r) in enumerate(
        results[(results['count']>1)&(results.Trait != 'Eczema')].iterrows()):
    results.loc[ind, 'gene'] += ''.join(['*']*(i+1))
results['Top TF (num)'] = ['{} ({:d})'.format(g, n)
        for g,n in zip(results.gene, results['count'])]
results['Sign'] = ['+' if r>0 else r'\---' for r in results.rf]

results.enrichment = ['{:.0%}'.format(x).replace('%',r'\%') for x in results.enrichment]
results.cell_line = [x + ' ()' for x in results.cell_line]

results = results[['Trait','Top TF (num)','rf','p','q',r'Min. \# sites']]

# print as tex
results.rename(
    columns={
        'p':'$p$',
        'q':'$q$',
        'rf':'$r_f$',
        }).to_latex(outfile,
                index=False,
                column_format='llrccc',
                bold_rows=True,
                escape=False)

# print as xlsx
results.to_excel(writer_xlsx, 'Table 1',
                index=False)
writer_xlsx.save()
