from __future__ import print_function, division
import pandas as pd
import os

sldp = '/groups/price/yakir/sldp'
annotinfodir = sldp + '/0.annotsummary/'

# phenotype names
phenotypes = {
        'ALZ' : 'Alzheimer\'s',
        'CD' : 'Crohn\'s',
        'UC' : 'Ulcerative colitis',
        'PASS_Anorexia' : 'Anorexia',
        'PASS_Autism' : 'Autism',
        'PASS_Celiac' : 'Celiac',
        'PASS_Coronary_Artery_Disease' : 'Coronary artery disease',
        'PASS_ENIGMA2_MeanPutamen' : 'Mean putamen volume',
        'PASS_Fasting_Glucose' : 'Fasting glucose',
        'PASS_HbA1C' : 'Hemoglobin A1C',
        'PASS_HDL' : 'HDL',
        'PASS_LDL' : 'LDL',
        'PASS_Lupus' : 'Lupus',
        'PASS_Rheumatoid_Arthritis' : 'Rheumatoid arthritis',
        'PASS_Schizophrenia' : 'Schizophrenia',
        'PASS_Primary_biliary_cirrhosis' : 'Primary biliary cirrhosis',
        'UKB_460K.blood_EOSINOPHIL_COUNT' : 'Eosinophil count',
        'UKB_460K.blood_PLATELET_COUNT' : 'Platelet count',
        'UKB_460K.blood_RBC_DISTRIB_WIDTH' : 'RBC distribution width',
        'UKB_460K.blood_RED_COUNT' : 'RBC count',
        'UKB_460K.blood_WHITE_COUNT' : 'WBC count',
        'UKB_460K.blood_MONOCYTE_COUNT' : 'Monocyte count',
        'UKB_460K.blood_LYMPHOCYTE_COUNT' : 'Lymphocyte count',
        'UKB_460K.blood_HIGH_LIGHT_SCATTER_RETICULOCYTE_COUNT' : 'Reticulocyte count',
        'UKB_460K.bmd_HEEL_TSCOREz' : 'Heel T score',
        'UKB_460K.body_BALDING1' : 'Balding',
        'UKB_460K.body_BMIz' : 'BMI',
        'UKB_460K.body_HEIGHTz' : 'Height',
        'UKB_460K.body_WHRadjBMIz' : 'Waist-hip ratio (BMI adj)',
        'UKB_460K.bp_SYSTOLICadjMEDz' : 'Systolic BP (med adj)',
        'UKB_460K.cov_EDU_YEARS' : 'Years of ed.',
        'UKB_460K.cov_SMOKING_STATUS' : 'Smoking status',
        'UKB_460K.disease_AID_SURE' : 'Auto-immune disease',
        'UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED' : 'Eczema',
        'UKB_460K.disease_DERMATOLOGY' : 'Dermatologic disease',
        'UKB_460K.disease_HI_CHOL_SELF_REP' : 'High cholesterol (self rep)',
        'UKB_460K.disease_HYPOTHYROIDISM_SELF_REP' : 'Hypothyroidism (self rep)',
        'UKB_460K.disease_RESPIRATORY_ENT' : 'Respiratory disease',
        'UKB_460K.disease_T2D' : 'Type 2 diabetes',
        'UKB_460K.lung_FEV1FVCzSMOKE' : 'FEV1/FVC (smoking adj)',
        'UKB_460K.lung_FVCzSMOKE' : 'FVC (smoking adj)',
        'UKB_460K.mental_NEUROTICISM' : 'Neuroticism',
        'UKB_460K.pigment_HAIR' : 'Hair color',
        'UKB_460K.pigment_SUNBURN' : 'Sunburn',
        'UKB_460K.repro_MENARCHE_AGE' : 'Age at menarche',
        'UKB_460K.repro_MENOPAUSE_AGE' : 'Age at menopause',

        'BP_mono_gene_nor_combat_peer_10' : 'GE (mono)',
        'BP_mono_K27AC_log2rpm_peer_10' : 'K27ac (mono)',
        'BP_mono_K4ME1_log2rpm_peer_10' : 'K4me1 (mono)',
        'BP_mono_meth_M_peer_10' : 'methylation (mono)',
        'BP_neut_gene_nor_combat_peer_10' : 'GE (neut)',
        'BP_neut_K27AC_log2rpm_peer_10' : 'K27ac (neut)',
        'BP_neut_K4ME1_log2rpm_peer_10' : 'K4me1 (neut)',
        'BP_neut_meth_M_peer_10' : 'methylation (neut)',
        'BP_tcel_gene_nor_combat_peer_10' : 'GE (T)',
        'BP_tcel_K27AC_log2rpm_peer_10' : 'K27ac (T)',
        'BP_tcel_K4ME1_log2rpm_peer_10' : 'K4me1 (T)',
        'BP_tcel_meth_M_peer_10' : 'methylation (T)',
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
