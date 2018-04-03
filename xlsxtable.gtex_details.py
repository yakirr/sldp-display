from __future__ import print_function, division
import pandas as pd
import numpy as np
import os
import ypy.fs as fs
from plot import params, results_overview, tables, info

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldprev+'/1.basset1tfs_p12/old_results/molecular_gtexv7_tissues/'
sumstats = params.sumstats + 'p12.molecular_gtexv7_tissues.phenos'
outname_list = me+'/out/xlsxtable.gtex_traitlist.xlsx'
writer_list = pd.ExcelWriter(outname_list)
outname_details = me+'/out/xlsxtable.gtex_results.xlsx'
writer_details = pd.ExcelWriter(outname_details)


## part 0
print('PART 0')
phenos = pd.read_csv(sumstats, header=None, names=['pheno']).set_index('pheno')
for p in phenos.index:
    print('getting info for', p)
    pinfo = pd.read_csv(params.sumstats + 'processed/' + p + '.KG3.95/info', sep='\t')
    phenos.loc[p, 'M'] = sum([
            pd.read_csv(params.sumstats + 'processed/' + p +'.KG3.95/' + str(c) + '.pss.gz',
                sep='\t').Winv_ahat_I.notnull().sum()
            for c in range(1,23)])
phenos['Trait'] = [info.phenotypes[p] for p in phenos.index]
phenos[['Trait', 'M']].to_excel(writer_list, 'Traits',
        index=False)

def format_and_print(passed, sheet_name):
    # format
    tables.format_sldp_results(passed, latex=False)
    passed.Trait = passed.Trait.str.replace('GE','gene expression')
    passed['activating'] = passed.uniprot_activator & ~passed.uniprot_repressor
    # print
    passed[['Trait','gene','cell_line','origin','rf','p',
        ]].rename(
        columns={
            'origin':'Lab',
            'cell_line':'Cell line',
            'gene':'TF',
            'rf':'$r_f$',
            'p':'$p$'}).to_excel(writer_details, sheet_name,
                    index=False)

## part a: unconditional GTEx results
results = results_overview.init(
        [indir+'all.gwresults'],
        [indir+'fdr5.gwresults'],
        'gene')
results['absrf'] = np.abs(results.rf)
results = results.sort_values(by='absrf', ascending=False)
passed = results[results.passed].copy()
format_and_print(passed, 'A. unconditional results')

## part b: conditional results
results = results_overview.init(
        [indir+'all.gwresults'],
        [indir+'fdr5_tissuespecific.gwresults'],
        'gene')
results['absrf'] = np.abs(results.rf)
results = results.sort_values(by='absrf', ascending=False)
passed = results[results.passed].copy()
format_and_print(passed, 'B. conditional results')


## save
writer_details.save()
writer_list.save()
