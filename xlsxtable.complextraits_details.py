from __future__ import print_function, division
import pandas as pd
import os
import pyutils.fs as fs
from plot import params, results_overview, tables, info

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/7.p9_a9/compiled_results/'
sumstats = params.sumstats + 'p9.complex.phenos'
outname_list = me+'/out/xlsxtable.complextraits_list.xlsx'
writer_list = pd.ExcelWriter(outname_list)
outname_details = me+'/out/xlsxtable.complextraits_details.xlsx'
writer_details = pd.ExcelWriter(outname_details)


## part 0
print('PART 0')
phenos = pd.read_csv(sumstats, header=None, names=['pheno']).set_index('pheno')
for p in phenos.index:
    print('getting info for', p)
    pinfo = pd.read_csv(params.sumstats + 'processed/' + p + '.KG3.95/info', sep='\t')
    phenos.loc[p, 'h2g'] = pinfo.h2g[0]
    phenos.loc[p, 'avg N'] = pinfo.Nbar[0]
    phenos.loc[p, 'M'] = sum([
            pd.read_csv(params.sumstats + 'processed/' + p +'.KG3.95/' + str(c) + '.pss.gz',
                sep='\t').Winv_ahat_I.notnull().sum()
            for c in range(1,23)])
phenos['Trait'] = [info.phenotypes[p] for p in phenos.index]
phenos[['Trait', 'avg N', 'M', 'h2g']].to_excel(writer_list, 'Traits',
        index=False)

## part a
print('PART A')
# get data
results = results_overview.init(
        [indir+'/p9.complex.all'],
        [indir+'/p9.complex.fdr5'],
        'gene').sort_values('sf_p')
passed = results[results.passed].copy()
# format
tables.format_sldp_results(passed, latex=False)
# print
passed[['Trait','gene','cell_line','origin','r_f','sf_p']].rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'gene':'TF',
        'r_f':'$r_f$',
        'sf_p':'$p$'}).to_excel(writer_details, 'A. TF associations',
                index=False)

## part b
print('PART B')
# get data
results = results_overview.init(
        [indir+'../../8.p9complex_a9minor/compiled_results/a9minor_maf5/p9.complex.all'],
        [indir+'../../8.p9complex_a9minor/compiled_results/a9minor_maf5/p9.complex.fdr5'],
        'gene').sort_values('sf_p')
passed = results[results.passed].copy()
# format
tables.format_sldp_results(passed, latex=False)
# print
passed[['Trait','gene','cell_line','origin','r_f','sf_p']].rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'gene':'TF',
        'r_f':'$r_f$',
        'sf_p':'$p$'}).to_excel(writer_details, 'B. Minor allele associations',
                index=False)

## part c
print('PART C')
results = results_overview.init(
        [indir+'/p9.complex.fdr5.indep'],
        [indir+'/p9.complex.fdr5.indep'],
        'cell_line').sort_values('sf_p')

# format
tables.format_sldp_results(results, latex=False)
results['TF (num)'] = ['{} ({:d})'.format(g, n)
        for g,n in zip(results.gene, results['count'])]
results.enrichment = ['{:.0%}'.format(x) for x in results.enrichment]
# print
results[['Trait','gene','cell_line','origin','r_f','sf_p','enrichment']].rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'gene':'TF',
        'sf_p':'p'}).to_excel(writer_details, 'C. unsigned enrichment',
                index=False)

## save
writer_list.save()
writer_details.save()
