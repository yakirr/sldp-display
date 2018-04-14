from __future__ import print_function, division
import pandas as pd
import os
import ypy.fs as fs
from plot import params, results_overview, tables, info

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldprev+'/1.basset1tfs_p12/old_results/complex/'
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
passed = results_overview.init(
        [indir+'/fdr5.sub_ng.gwresults'],
        [indir+'/fdr5.sub_ng.gwresults'],
        'gene').sort_values('p')
# format
tables.format_sldp_results(passed, latex=False)
# print
passed[['Trait','gene','cell_line','origin','rf','p','q']].rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'gene':'TF',
        'rf':'$rf$',
        'p':'$p$',
        'q':'$q$'}).to_excel(writer_details, 'A. TF associations',
                index=False)

## part b
print('PART B')
# get data
passed = results_overview.init(
        [indir+'/minor.fdr5.sub_ng.gwresults'],
        [indir+'/minor.fdr5.sub_ng.gwresults'],
        'gene').sort_values('p')
# format
tables.format_sldp_results(passed, latex=False)
# print
passed[['Trait','gene','cell_line','origin','rf','p']].rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'gene':'TF',
        'rf':'$rf$',
        'p':'$p$'}).to_excel(writer_details, 'B. Minor allele associations',
                index=False)

## part c
print('PART C')
results = results_overview.init(
        [indir+'/fdr5_indep.sub_ng.gwresults'],
        [indir+'/fdr5_indep.sub_ng.gwresults'],
        'cell_line').sort_values('p')

# format
tables.format_sldp_results(results, latex=False)
results['TF (num)'] = ['{} ({:d})'.format(g, n)
        for g,n in zip(results.gene, results['count'])]
results.enrichment = ['{:.0%}'.format(x) for x in results.enrichment]
# print
results[['Trait','gene','cell_line','origin','rf','p','enrichment']].rename(
    columns={
        'origin':'Lab',
        'cell_line':'Cell line',
        'gene':'TF',
        'p':'p'}).to_excel(writer_details, 'C. heritability explained',
                index=False)

## save
writer_list.save()
writer_details.save()
