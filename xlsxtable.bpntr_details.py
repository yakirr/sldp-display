from __future__ import print_function, division
import pandas as pd
import numpy as np
import os
import ypy.fs as fs
from plot import params, results_overview, tables, info

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldprev+'/1.basset1tfs_p12/'
sumstats = params.sumstats + 'p12.molecular_bpntr.phenos'
outname_list = me+'/out/xlsxtable.bpntr_traitlist.xlsx'
writer_list = pd.ExcelWriter(outname_list)
outname_details = me+'/out/xlsxtable.bpntr_results.xlsx'
writer_details = pd.ExcelWriter(outname_details)


## part 0
# print('PART 0')
# phenos = pd.read_csv(sumstats, header=None, names=['pheno']).set_index('pheno')
# for p in phenos.index:
#     print('getting info for', p)
#     pinfo = pd.read_csv(params.sumstats + 'processed/' + p + '.KG3.95/info', sep='\t')
#     phenos.loc[p, 'M'] = sum([
#             pd.read_csv(params.sumstats + 'processed/' + p +'.KG3.95/' + str(c) + '.pss.gz',
#                 sep='\t').Winv_ahat_I.notnull().sum()
#             for c in range(1,23)])
# phenos['Trait'] = [info.phenotypes[p] for p in phenos.index]
# phenos.Trait = phenos.Trait.str.replace('GE','gene expression')
# phenos[['Trait', 'M']].to_excel(writer_list, 'Traits',
#         index=False)

def format_and_print(passed, sheet_name):
    # format
    tables.format_sldp_results(passed, latex=False)
    passed.Trait = passed.Trait.str.replace('GE','gene expression')
    passed['activating'] = passed.uniprot_activator & ~passed.uniprot_repressor
    # print
    passed[['Trait','gene','cell_line','origin','rf','p',
        'activating','repressing','ambiguous']].rename(
        columns={
            'origin':'Lab',
            'cell_line':'Cell line',
            'gene':'TF',
            'rf':'$r_f$',
            'p':'$p$'}).to_excel(writer_details, sheet_name,
                    index=False)

## part a: blueprint GE
results = results_overview.init(
        [indir+'molecular_BP/all.gwresults'],
        [indir+'molecular_BP/fdr5.gwresults'],
        'gene')
results['absrf'] = np.abs(results.rf)
results = results.sort_values(by='absrf', ascending=False)
results['activating'] = results.uniprot_activator & ~results.uniprot_repressor
results['repressing'] = ~results.uniprot_activator & results.uniprot_repressor
results['ambiguous'] = ~results.activating & ~results.repressing
passed = results[results.passed & results.pheno.str.contains('gene_')].copy()
format_and_print(passed, 'A. gene expression (BLUEPRINT)')
print(len(passed[['Trait','gene']].drop_duplicates()), 'distinct trait gene pairs')
print()

## part b: ntr GE
results = results_overview.init(
        [indir+'/molecular_NTR/all.gwresults'],
        [indir+'/molecular_NTR/fdr5.gwresults'],
        'gene')
results['absrf'] = np.abs(results.rf)
results = results.sort_values(by='absrf', ascending=False)
results['activating'] = results.uniprot_activator & ~results.uniprot_repressor
results['repressing'] = ~results.uniprot_activator & results.uniprot_repressor
results['ambiguous'] = ~results.activating & ~results.repressing
passed = results[results.passed & results.pheno.str.contains('NTR')].copy()
format_and_print(passed, 'B. gene expression (NTR)')
print(len(passed[['Trait','gene']].drop_duplicates()), 'distinct trait gene pairs')
print()

## part c
results = results_overview.init(
        [indir+'molecular_BP/all.gwresults', indir+'/molecular_NTR/all.gwresults'],
        [indir+'molecular_BP/fdr5.gwresults', indir+'/molecular_NTR/fdr5.gwresults'],
        'gene')
results['absrf'] = np.abs(results.rf)
results = results.sort_values(by='absrf', ascending=False)
results['activating'] = results.uniprot_activator & ~results.uniprot_repressor
results['repressing'] = ~results.uniprot_activator & results.uniprot_repressor
results['ambiguous'] = ~results.activating & ~results.repressing
x = results[results.pheno.str.contains('neut') &
        results.pheno.str.contains('gene')].sort_values('annot')
y = results[results.pheno.str.contains('NTR')].sort_values('annot')
both = x[['origin','cell_line','gene','activating','repressing','ambiguous','z']]
both.loc[:,'z_NTR'] = y.z.values
both[['gene','cell_line','origin','z','z_NTR','activating','repressing','ambiguous']].rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'gene':'TF',
        'z':'z-score (BLUEPRINT)',
        'z_NTR':'z-score (NTR)'}).to_excel(writer_details, 'C. BLUEPRINT vs NTR',
                index=False)


# part d: bp k4me1
results = results_overview.init(
        [indir+'molecular_BP/all.gwresults'],
        [indir+'molecular_BP/fdr5.gwresults'],
        'gene')
results['absrf'] = np.abs(results.rf)
results = results.sort_values(by='absrf', ascending=False)
results['activating'] = results.uniprot_activator & ~results.uniprot_repressor
results['repressing'] = ~results.uniprot_activator & results.uniprot_repressor
results['ambiguous'] = ~results.activating & ~results.repressing
passed = results[results.passed & results.pheno.str.contains('K4ME1')].copy()
format_and_print(passed, 'D. K4me1 (BLUEPRINT)')
print(len(passed[['Trait','gene']].drop_duplicates()), 'distinct trait gene pairs')
print()

# part e: bp k27ac
results = results_overview.init(
        [indir+'molecular_BP/all.gwresults'],
        [indir+'molecular_BP/fdr5.gwresults'],
        'gene')
results['absrf'] = np.abs(results.rf)
results = results.sort_values(by='absrf', ascending=False)
results['activating'] = results.uniprot_activator & ~results.uniprot_repressor
results['repressing'] = ~results.uniprot_activator & results.uniprot_repressor
results['ambiguous'] = ~results.activating & ~results.repressing
passed = results[results.passed & results.pheno.str.contains('K27AC')].copy()
format_and_print(passed, 'E. K27ac (BLUEPRINT)')
print(len(passed[['Trait','gene']].drop_duplicates()), 'distinct trait gene pairs')
print()


## save
# writer_details.save()
# writer_list.save()
