from __future__ import print_function, division
import pandas as pd
import os
import ypy.fs as fs
from plot import params, results_overview, tables, info

me = os.path.dirname(os.path.abspath(__file__))
enrichments_file = params.sldprev+'/6.complex_msigdb/fdr5.enrichments'
outname_numerical = me+'/out/xlsxtable.enrichment_details.xlsx'
writer_numerical = pd.ExcelWriter(outname_numerical)

# read in enrichment results and split into traits and annotations and format p-values
enrichments = pd.read_csv(enrichments_file, sep='\t')
enrichments['Trait'] = ['.'.join(t.split('.')[:-1]) for t in enrichments.target]
enrichments['Annotation'] = enrichments.target.str.split('.').str.get(-1)

# sort enrichments
enrichments['xabs'] = enrichments.mean_in_wgt/enrichments.mean_out_wgt
sortkey=['Trait','Annotation','q','xabs']
sortkey=['q','xabs']
print('sorting by', sortkey)
enrichments = enrichments.sort_values(by=sortkey, ascending=[True,False])

# format limiting p-values
enrichments.loc[enrichments.p_wgt == 0, 'p_wgt'] = '<1e-6'
# enrichments.loc[enrichments.q == 0, 'q'] = '-'
enrichments.loc[enrichments.q == 0, 'q'] = [ 1e-6*nt for nt in
        enrichments.loc[enrichments.q == 0, 'numtests']]

# add annotation metadata
enrichments = pd.merge(enrichments, info.annots[['annot','origin','cell_line','gene']],
        left_on='Annotation', right_on='annot',
        how='left')

enrichments = enrichments[['Trait','gene','cell_line','origin',
                            'geneset','mean_in_wgt','mean_out_wgt',
                            'std_in_wgt','std_out_wgt','p_wgt','q']].rename(
            columns={
                'gene':'TF',
                'cell_line':'Cell line',
                'origin':'Lab',
                'geneset':'Gene set',
                'mean_in_wgt':'Avg covariance in set',
                'mean_out_wgt':'Avg covariance outside set',
                'std_in_wgt':'Std error(Avg in set)',
                'std_out_wgt':'Std error(Avg outside set)',
                'p_wgt':'p',
                })

enrichments.to_excel(
        writer_numerical, 'A. MSigDB enrichments',
        index=False)

## save
writer_numerical.save()
