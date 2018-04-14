from __future__ import print_function, division
import pandas as pd
import numpy as np
import os
import ypy.fs as fs
from plot import params, results_overview, tables, info

me = os.path.dirname(os.path.abspath(__file__))

def format_and_print(passed, sheet_name, writer):
    # format
    tables.format_sldp_results(passed, latex=False)
    if len(passed) > 0:
        passed.Trait = passed.Trait.str.replace('GE','gene expression')
    # print
    passed[['Trait','gene','cell_line','origin','rf','p',
        ]].rename(
        columns={
            'origin':'Lab',
            'cell_line':'Cell line',
            'gene':'TF',
            'rf':'$r_f$',
            'p':'$p$'}).to_excel(writer, sheet_name,
                    index=False)

def format_and_print_gtrd(passed, sheet_name, writer):
    # format
    tables.format_sldp_results(passed, latex=False)
    if len(passed) > 0:
        passed.Trait = passed.Trait.str.replace('GE','gene expression')
    # print
    passed[['Trait','GTRD annotation','rf','p',
        ]].rename(
        columns={
            'rf':'$r_f$',
            'p':'$p$'}).to_excel(writer, sheet_name,
                    index=False)

traitsets = [
        ('molecular_BP', 'A. BLUEPRINT'),
        ('molecular_NTR', 'B. NTR'),
        ('molecular_gtexv7_tissues', 'C. GTEx'),
        ('complex', 'D. Complex traits'),
        ]

## deepsea
indir = params.sldprev + '/3.deepsea_p12/basset1qc/'
outname = me+'/out/xlsxtable.deepsea_results.xlsx'
writer = pd.ExcelWriter(outname)
for traits, sheet_name in traitsets:
    results = results_overview.init(
            [indir+traits+'/fdr5.gwresults'],
            [indir+traits+'/fdr5.gwresults'],
            'gene')
    results['absrf'] = np.abs(results.rf)
    results = results.sort_values(by='absrf', ascending=False)
    format_and_print(results, sheet_name, writer)
writer.save()

## hocomoco
indir = params.sldprev + '/2.hocomotifwbasset_p12/expdiff_sum/'
outname = me+'/out/xlsxtable.motif_results.xlsx'
writer = pd.ExcelWriter(outname)
for traits, sheet_name in traitsets:
    results = results_overview.init(
            [indir+traits+'/fdr5.gwresults'],
            [indir+traits+'/fdr5.gwresults'],
            'gene')
    results['absrf'] = np.abs(results.rf)
    results = results.sort_values(by='absrf', ascending=False)
    format_and_print(results, sheet_name, writer)
writer.save()

## gtrd
indir = params.sldprev + '/5.gtrdbasset2tfs_p12/gtrd/'
outname = me+'/out/xlsxtable.gtrd_results.xlsx'
writer = pd.ExcelWriter(outname)
for traits, sheet_name in traitsets:
    results = pd.read_csv(indir+traits+'/fdr5.gwresults', sep='\t')
    results['Trait'] = [info.phenotypes[p] for p in results.pheno]
    results['GTRD annotation'] = [a.replace('_GTRD.R','') for a in results.annot]
    results['absrf'] = np.abs(results.rf)
    results = results.sort_values(by='absrf', ascending=False)
    format_and_print_gtrd(results, sheet_name, writer)
writer.save()
