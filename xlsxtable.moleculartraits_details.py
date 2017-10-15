from __future__ import print_function, division
import pandas as pd
import os
import pyutils.fs as fs
from plot import params, results_overview, tables, info

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/7.p9_a9/compiled_results/'
sumstats = params.sumstats + 'p9.molecular.phenos'
outname_list = me+'/out/xlsxtable.moleculartraits_list.xlsx'
writer_list = pd.ExcelWriter(outname_list)
outname_details = me+'/out/xlsxtable.moleculartraits_details.xlsx'
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
phenos.Trait = phenos.Trait.str.replace('GE','gene expression')
phenos[['Trait', 'M']].to_excel(writer_list, 'Traits',
        index=False)

def format_and_print(passed, sheet_name):
    # format
    tables.format_sldp_results(passed, latex=False)
    passed.Trait = passed.Trait.str.replace('GE','gene expression')
    passed['activating'] = passed.uniprot_activator & ~passed.uniprot_repressor
    # print
    passed[['Trait','gene','cell_line','origin','r_f','sf_p','activating']].rename(
        columns={
            'origin':'Lab',
            'cell_line':'Cell line',
            'gene':'TF',
            'r_f':'$r_f$',
            'sf_p':'$p$'}).to_excel(writer_details, sheet_name,
                    index=False)

## part a
results = results_overview.init(
        [indir+'/p9.molecular.all'],
        [indir+'/p9.molecular.fdr5'],
        'gene').sort_values('sf_p')
passed = results[results.passed & results.pheno.str.contains('gene_')].copy()
format_and_print(passed, 'A. gene expression (BLUEPRINT)')

## part b
results = results_overview.init(
        [indir+'/p9.molecular.all'],
        [indir+'/p9.molecular.fdr5'],
        'gene').sort_values('sf_p')
passed = results[results.passed & results.pheno.str.contains('NTR')].copy()
format_and_print(passed, 'B. gene expression (NTR)')

## part c
results = results_overview.init(
        [indir+'/p9.molecular.all'],
        [indir+'/p9.molecular.fdr5'],
        'gene').sort_values('sf_p')
results['activating'] = results.uniprot_activator & ~results.uniprot_repressor
x = results[results.pheno.str.contains('neut') &
        results.pheno.str.contains('gene')].sort_values('annot')
y = results[results.pheno.str.contains('NTR')].sort_values('annot')
both = x[['origin','cell_line','gene','activating','sf_z']]
both.loc[:,'sf_z_NTR'] = y.sf_z.values
both[['gene','cell_line','origin','sf_z','sf_z_NTR','activating']].rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'gene':'TF',
        'sf_z':'z-score (BLUEPRINT)',
        'sf_z_NTR':'z-score (NTR)'}).to_excel(writer_details, 'C. BLUEPRINT vs NTR',
                index=False)


## part d
results = results_overview.init(
        [indir+'/p9.molecular.all'],
        [indir+'/p9.molecular.fdr5'],
        'gene').sort_values('sf_p')
passed = results[results.passed & results.pheno.str.contains('K4ME1')].copy()
format_and_print(passed, 'D. K4me1 (BLUEPRINT)')

# part e
results = results_overview.init(
        [indir+'/p9.molecular.all'],
        [indir+'/p9.molecular.fdr5'],
        'gene').sort_values('sf_p')
passed = results[results.passed & results.pheno.str.contains('K27AC')].copy()
format_and_print(passed, 'E. K27ac (BLUEPRINT)')

# part f
results = results_overview.init(
        [indir+'/p9.molecular.all'],
        [indir+'/p9.molecular.fdr10'],
        'gene').sort_values('sf_p')
passed = results[results.passed & results.pheno.str.contains('meth')].copy()
format_and_print(passed, 'F. methyl (BLUEPRINT; FDR <10%)')



## save
writer_details.save()
writer_list.save()
