from __future__ import print_function, division
import pandas as pd
import os

sldp = '../sldp.sub_ng/'
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

        'gtexv7h2g50sqrt_Adipose_Subcutaneous' : 'Adipose (subcut.)',
        'gtexv7h2g50sqrt_Adipose_Visceral_Omentum' : 'Adipose (visc. om.)',
        'gtexv7h2g50sqrt_Adrenal_Gland' : 'Adrenal gl.',
        'gtexv7h2g50sqrt_Artery_Aorta' : 'SM Artery (aorta)',
        'gtexv7h2g50sqrt_Artery_Coronary' : 'SM Artery (coron.)',
        'gtexv7h2g50sqrt_Artery_Tibial' : 'SM Artery (tibial)',
        'gtexv7h2g50sqrt_Brain_Amygdala' : 'Brain (amygdala)',
        'gtexv7h2g50sqrt_Brain_Anterior_cingulate_cortex_BA24' : 'Brain (ant. cing. cort.)',
        'gtexv7h2g50sqrt_Brain_Caudate_basal_ganglia' : 'Brain (caudate)',
        'gtexv7h2g50sqrt_Brain_Cerebellar_Hemisphere' : 'Brain (cerebel. hem.)',
        'gtexv7h2g50sqrt_Brain_Cerebellum' : 'Brain (cerebel.)',
        'gtexv7h2g50sqrt_Brain_Cortex' : 'Brain (cort.)',
        'gtexv7h2g50sqrt_Brain_Frontal_Cortex_BA9' : 'Brain (frontal cort.)',
        'gtexv7h2g50sqrt_Brain_Hippocampus' : 'Brain (hipp.)',
        'gtexv7h2g50sqrt_Brain_Hypothalamus' : 'Brain (hypothal.)',
        'gtexv7h2g50sqrt_Brain_Nucleus_accumbens_basal_ganglia' : 'Brain (nuc. acum.)',
        'gtexv7h2g50sqrt_Brain_Putamen_basal_ganglia' : 'Brain (putamen)',
        'gtexv7h2g50sqrt_Brain_Spinal_cord_cervical_c-1' : 'Brain (C1)',
        'gtexv7h2g50sqrt_Brain_Substantia_nigra' : 'Brain (sub. nig.)',
        'gtexv7h2g50sqrt_Breast_Mammary_Tissue' : 'Breast (mammary)',
        'gtexv7h2g50sqrt_Cells_EBV-transformed_lymphocytes' : 'Lymphocytes',
        'gtexv7h2g50sqrt_Cells_Transformed_fibroblasts' : 'Fibroblasts',
        'gtexv7h2g50sqrt_Colon_Sigmoid' : 'GI (sig. col.)',
        'gtexv7h2g50sqrt_Colon_Transverse' : 'GI (trans. col.)',
        'gtexv7h2g50sqrt_Esophagus_Gastroesophageal_Junction' : 'GI (gastroesoph. j.)',
        'gtexv7h2g50sqrt_Esophagus_Mucosa' : 'GI (esophagus muc.)',
        'gtexv7h2g50sqrt_Esophagus_Muscularis' : 'GI (esophagus musc.)',
        'gtexv7h2g50sqrt_Heart_Atrial_Appendage' : 'SM heart (atrial appen.)',
        'gtexv7h2g50sqrt_Heart_Left_Ventricle' : 'SM heart (l. vent.)',
        'gtexv7h2g50sqrt_Liver' : 'GI (liver)',
        'gtexv7h2g50sqrt_Lung' : 'Lung',
        'gtexv7h2g50sqrt_Minor_Salivary_Gland' : 'Min. salivary gl.',
        'gtexv7h2g50sqrt_Muscle_Skeletal' : 'Skeletal muscle',
        'gtexv7h2g50sqrt_Nerve_Tibial' : 'Tibial nerve',
        'gtexv7h2g50sqrt_Ovary' : 'Ovary',
        'gtexv7h2g50sqrt_Pancreas' : 'Pancreas',
        'gtexv7h2g50sqrt_Pituitary' : 'Pituitary',
        'gtexv7h2g50sqrt_Prostate' : 'Prostate',
        'gtexv7h2g50sqrt_Skin_Not_Sun_Exposed_Suprapubic' : 'Skin (not sun exp.)',
        'gtexv7h2g50sqrt_Skin_Sun_Exposed_Lower_leg' : 'Skin (sun exp.)',
        'gtexv7h2g50sqrt_Small_Intestine_Terminal_Ileum' : 'GI (sm. int.)',
        'gtexv7h2g50sqrt_Spleen' : 'Spleen',
        'gtexv7h2g50sqrt_Stomach' : 'GI (stomach)',
        'gtexv7h2g50sqrt_Testis' : 'Testis',
        'gtexv7h2g50sqrt_Thyroid' : 'Thyroid',
        'gtexv7h2g50sqrt_Uterus' : 'Uterus',
        'gtexv7h2g50sqrt_Vagina' : 'Vagina',
        'gtexv7h2g50sqrt_Whole_Blood' : 'Whole blood',

        'BPh2g50sqrt_mono_gene_nor_combat' : 'GE (mono)',
        'BPh2g50sqrt_mono_K27AC_log2rpm' : 'K27ac (mono)',
        'BPh2g50sqrt_mono_K4ME1_log2rpm' : 'K4me1 (mono)',
        'BPh2g50sqrt_mono_meth_M' : 'methylation (mono)',
        'BPh2g50sqrt_neut_gene_nor_combat' : 'GE (neut)',
        'BPh2g50sqrt_neut_K27AC_log2rpm' : 'K27ac (neut)',
        'BPh2g50sqrt_neut_K4ME1_log2rpm' : 'K4me1 (neut)',
        'BPh2g50sqrt_neut_meth_M' : 'methylation (neut)',
        'BPh2g50sqrt_tcel_gene_nor_combat' : 'GE (T)',
        'BPh2g50sqrt_tcel_K27AC_log2rpm' : 'K27ac (T)',
        'BPh2g50sqrt_tcel_K4ME1_log2rpm' : 'K4me1 (T)',
        'BPh2g50sqrt_tcel_meth_M' : 'methylation (T)',
        'NTRh2g50sqrt_blood_gene' : 'GE (NTR)'
        }

# compile metadata about experiments
experiment_to_gene = pd.read_csv(annotinfodir + 'experiment_to_gene.tsv', sep='\t')
gene_to_uniprotinfo = pd.read_csv(annotinfodir + 'gene_to_uniprotinfo.tsv', sep='\t')
gene_to_category = pd.read_csv(annotinfodir + 'gene_to_category.tsv', sep=',')
gene_to_ensembl = pd.read_csv(annotinfodir + 'gene_to_ensgid.txt', delim_whitespace=True)
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
