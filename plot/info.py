from __future__ import print_function, division
import pandas as pd
import os

sldp = '/groups/price/yakir/sldp'

# phenotype names
phenotypes = {
        'CD' : 'CD',
        'IBD' : 'IBD',
        'PASS_Anorexia' : 'Anorexia',
        'PASS_HDL' : 'HDL',
        'PASS_Lupus' : 'Lupus',
        'UKBiobank_Eczema' : 'Eczema',
        'BP_mono_gene_nor_combat_peer_10' : 'GE (mono)',
        'BP_mono_K4ME1_log2rpm_peer_10' : 'K4me1 (mono)',
        'BP_neut_gene_nor_combat_peer_10' : 'GE (neut)',
        'BP_neut_K27AC_log2rpm_peer_10' : 'K27ac (neut)',
        'BP_neut_K4ME1_log2rpm_peer_10' : 'K4me1 (neut)',
        'geneexp_total_NTR' : 'GE (NTR)'
        }

# annotations with metadata
annots = pd.read_csv(sldp + '/0.annotsummary/annotsummary.tsv', sep='\t').rename(
        columns={'id':'annot'})
annots = annots[annots.supp >= 5000]

annots['bigexp'] = annots.experiment
annots.loc[annots.experiment.str.contains('CTCF') |
        annots.experiment.str.contains('RAD21'), 'bigexp'] = 'CTCF/RAD21'
annots.loc[annots.experiment.str.contains('POL2') |
        annots.experiment.str.contains('TBP') |
        annots.experiment.str.contains('TAF1'), 'bigexp'] = 'POL2/TBP/TAF1'
annots['tissue'] = 'non-blood'
annots.loc[annots.cell_line.str.contains('GM') |
        annots.cell_line.str.contains('K562') |
        annots.cell_line.str.contains('PBDE') |
        annots.cell_line.str.contains('HUVEC') |
        annots.cell_line.str.contains('RAJI'), 'tissue'] = 'blood'
