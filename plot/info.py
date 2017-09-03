from __future__ import print_function, division
import pandas as pd
import os

sldp = '/groups/price/yakir/sldp'
annotinfodir = sldp + '/0.annotsummary/'

# phenotype names
phenotypes = {
        'CD' : 'Crohn\'s',
        'PASS_Anorexia' : 'Anorexia',
        'PASS_HDL' : 'HDL',
        'PASS_Lupus' : 'Lupus',
        'UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED' : 'Eczema',
        'UKB_460K.cov_EDU_YEARS' : 'Years of Ed.',
        'BP_mono_gene_nor_combat_peer_10' : 'GE (mono)',
        'BP_mono_K27AC_log2rpm_peer_10' : 'K27ac (mono)',
        'BP_mono_K4ME1_log2rpm_peer_10' : 'K4me1 (mono)',
        'BP_neut_gene_nor_combat_peer_10' : 'GE (neut)',
        'BP_neut_K27AC_log2rpm_peer_10' : 'K27ac (neut)',
        'BP_neut_K4ME1_log2rpm_peer_10' : 'K4me1 (neut)',
        'BP_tcel_K4ME1_log2rpm_peer_10' : 'K4me1 (T)',
        'geneexp_total_NTR' : 'GE (NTR)'
        }

# compile metadata about experiments
experiment_to_gene = pd.read_csv(annotinfodir + 'experiment_to_gene.tsv', sep='\t')
gene_to_uniprotinfo = pd.read_csv(annotinfodir + 'gene_to_uniprotinfo.tsv', sep='\t')
gene_to_category = pd.read_csv(annotinfodir + 'gene_to_category.tsv', sep=',')
experiment_meta = pd.merge(experiment_to_gene,
        gene_to_uniprotinfo,
        on='gene', how='left')
experiment_meta = pd.merge(experiment_meta,
        gene_to_category,
        on='gene', how='left')
with open(annotinfodir + 'category_order.txt') as f:
    category_order = f.read().splitlines()

# read annotations with metadata and merge in experiment info
annots = pd.read_csv(annotinfodir + 'annotsummary.tsv', sep='\t').rename(
        columns={'id':'annot'})
annots = annots[annots.supp >= 5000]

annots = pd.merge(annots, experiment_meta, on='experiment', how='left')

# create description and gene group columns for annotations
annots['desc'] = annots.function
annots.loc[annots.function.isnull(), 'desc'] = \
        annots.loc[annots.function.isnull(), 'category']
annots.desc = pd.Categorical(annots.desc,
        categories=category_order,
        ordered=True)

annots['genegroup'] = annots.gene
annots.loc[annots.gene.str.contains('CTCF') |
        annots.gene.str.contains('RAD21'), 'genegroup'] = 'CTCF/RAD21'
annots.loc[annots.gene.str.contains('POL2') |
        annots.gene.str.contains('TBP') |
        annots.gene.str.contains('TAF1'), 'genegroup'] = 'POL2/TBP/TAF1'

# create tissue column for annotations
annots['tissue'] = 'non-blood'
annots.loc[annots.cell_line.str.contains('GM') |
        annots.cell_line.str.contains('K562') |
        annots.cell_line.str.contains('PBDE') |
        annots.cell_line.str.contains('HUVEC') |
        annots.cell_line.str.contains('RAJI'), 'tissue'] = 'blood'
