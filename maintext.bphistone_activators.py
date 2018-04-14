from __future__ import print_function, division
import plot.results_overview as ro
import scipy.stats as st

# for positive histone results
results = ro.init(['../sldp.rev1_ng/1.basset1tfs_p12/molecular_BP/all.gwresults'],
        ['../sldp.rev1_ng/1.basset1tfs_p12/molecular_BP/fdr5.gwresults'], 'desc')
results['sa'] = \
        results.uniprot_activator & ~results.uniprot_repressor
check = results[
        (results.pheno.str.contains('K27')|results.pheno.str.contains('K4'))&
        (results.rf>0) & results.passed]
print(len(check), check.sa.sum())

check = check.drop_duplicates(subset=['gene','pheno'])
n = check.sa.sum()
N = len(check)
p = results.sa.sum()/len(results)
print(N, p, n)
print(st.binom.sf(n-1, N, p))

