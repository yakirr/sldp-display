from __future__ import print_function, division
import pandas as pd
import numpy as np
import os
import ypy.fs as fs
from plot import params, results_overview, tables, info

me = os.path.dirname(os.path.abspath(__file__))
enrichments_file = params.sldprev+'/6.complex_msigdb/fdr5.enrichments'
corrected_enrichments_file = params.sldprev+'/7.complex_msigdb_repeat/fdr5.enrichments'
outname_numerical = me+'/out/xlsxtable.enrichment_details.compare_to_correction.raw.xlsx'
writer_numerical = pd.ExcelWriter(outname_numerical)

# read in original enrichment results and split into traits and annotations and format p-values
enrichments = pd.read_csv(enrichments_file, sep='\t')
enrichments['Trait'] = ['.'.join(t.split('.')[:-1]) for t in enrichments.target]
enrichments['Annotation'] = enrichments.target.str.split('.').str.get(-1)

# read in corrected results
corrected_enrichments = pd.read_csv(corrected_enrichments_file, sep='\t')
corrected_enrichments['Trait'] = ['.'.join(t.split('.')[:-1]) for t in corrected_enrichments.target]
corrected_enrichments['Annotation'] = corrected_enrichments.target.str.split('.').str.get(-1)

# join new p-values and q-values
enrichments = pd.merge(enrichments,
        corrected_enrichments[['target','geneset','p_wgt','q']].rename(
            {'p_wgt':'p_wgt_corrected','q':'q_corrected'}, axis='columns'),
        how='left',
        on=['target','geneset'])

# sort enrichments
enrichments['x'] = enrichments.mean_in_wgt/enrichments.mean_out_wgt
sortkey=['Trait','Annotation','q','x']
sortkey=['q','x']
print('sorting by', sortkey)
enrichments = enrichments.sort_values(by=sortkey, ascending=[True,False])

# format limiting p-values
enrichments.loc[enrichments.p_wgt == 0, 'p_wgt'] = '<1e-6'
enrichments.loc[enrichments.q == 0, 'q'] = '-'
enrichments.loc[enrichments.p_wgt_corrected == 0, 'p_wgt_corrected'] = '<1e-6'
enrichments.loc[np.isnan(enrichments.q_corrected), 'q_corrected'] = '>5%'
enrichments.loc[enrichments.q_corrected == 0, 'q_corrected'] = '-'

# add annotation metadata
enrichments = pd.merge(enrichments, info.annots[['annot','origin','cell_line','gene']],
        left_on='Annotation', right_on='annot',
        how='left')

enrichments = enrichments[['Trait','gene','cell_line','origin',
                            'geneset','x','mean_in_wgt','mean_out_wgt',
                            'std_in_wgt','std_out_wgt','p_wgt','q',
                            'p_wgt_corrected','q_corrected']].rename(
            columns={
                'gene':'TF',
                'cell_line':'Cell line',
                'origin':'Lab',
                'geneset':'Gene set',
                'x':'Enrichment',
                'mean_in_wgt':'Avg covariance in set',
                'mean_out_wgt':'Avg covariance outside set',
                'std_in_wgt':'Std error(Avg in set)',
                'std_out_wgt':'Std error(Avg outside set)',
                'p_wgt':'p',
                'p_wgt_corrected':'p_corrected',
                'q_corrected':'q_corrected',
                })

enrichments.to_excel(
        writer_numerical, 'B. All MSigDB enrichments',
        index=False)

## save
writer_numerical.save()
