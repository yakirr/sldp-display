from __future__ import print_function, division
from plot import params, info, results_overview

def format_p(p, latex=True):
    temp = '{:.1e}'.format(p)
    base = temp[:temp.index('e')]
    exp = str(int(temp[temp.index('e')+1:]))
    if latex:
        return '$'+base+'\\times 10^{'+exp+'}$'
    else:
        return base+'e'+exp

def format_sldp_results(results, latex=True):
    results['Trait'] = [info.phenotypes[p] for p in results.pheno]
    results.p = [format_p(x, latex) for x in results.p]
    results.rf = ['{:.1%}'.format(x).replace(
        '%',
        r'\%' if latex else r'%'
        ) for x in results.rf]
