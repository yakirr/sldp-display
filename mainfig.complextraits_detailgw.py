from __future__ import print_function, division
import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyutils.fs as fs
import statutils.vis as vis
from plot import params, results_overview

me = os.path.dirname(os.path.abspath(__file__))
indir = params.sldp+'/7.p9_a9/verboseresults/'
outname = me+'/out/mainfig.complextraits_detailgw.raw.pdf'

# set parameters
toplot = [
        # ('CD', 'SydhGm18951Pol2Iggmus', 'POL2'),
        # ('UKB_460K.disease_ALLERGY_ECZEMA_DIAGNOSED','UtaMcf7CtcfSerumstim', 'CTCF'),
        ('UKB_460K.cov_EDU_YEARS','HaibGm12878Bcl11aPcr1x', 'EDU', 'BCL11A', 'LCL',
            2, 60, 61, None),
        ('PASS_Lupus','SydhK562CtcfbIggrab', 'SLE', 'CTCF', 'K562', None, None, None, None),
        ('CD', 'SydhK562Irf1Ifng6h', 'CD', 'IRF1', 'K562', 5, 131.3, 133.49,
            '/groups/price/yakir/dpos/WEIGHTS/NTR.BLOOD.RNAARR/NTR.IRF1.wgt.RDat')
        ]
nrows = 3; ncols = 3

# set up figure
fig = plt.figure(figsize=(4,4))
gs = gridspec.GridSpec(nrows,ncols)

# make figure
for i,(pheno, annot, prettypheno, tf, cell_line, c, start, end, twas) in enumerate(toplot):
    print(pheno, tf)
    annotname = annot+',diff,'+annot+',-1.R'
    print(indir+pheno+'.'+annotname+'.toplot.gz')
    data = pd.read_csv(indir+pheno+'.'+annotname+'.toplot.gz', sep='\t')
    chunks = pd.read_csv(indir+pheno+'.'+annotname+'.chunks', sep='\t')
    mu = chunks.q.sum()/chunks.r.sum()

    # make plot of ahat vs Rv
    ax = plt.subplot(gs[i,0])
    ax.locator_params(nbins=4)
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.yaxis.offsetText.set_fontsize(params.labelfontsize)
    xs, ys = vis.scatter_b(
            data.Rv_resid, data.ahat_resid, binsize=4000, extreme_only=5, errorbars=True,
            ax=ax,
            color='purple',
            fmt='-o',
            linewidth=0.5,
            markersize=2)
    ax.plot([min(xs), max(xs)], [mu*min(xs), mu*max(xs)],
        c='gray', linewidth=0.8, dashes=[2,2])
    ax.tick_params(**params.tickprops)
    ax.set_xlabel(r'avg $Rv$ ('+tf+', '+cell_line+')', fontsize=params.labelfontsize)
    ax.set_ylabel(r'avg $\hat\alpha$ ('+prettypheno+')', fontsize=params.labelfontsize)

    # make plot of q
    ax = plt.subplot(gs[i,1])

    qsort = np.sort(chunks.q)
    qs = pd.DataFrame()
    qs['q'] = qsort
    half = np.argmin(qsort**2)
    percent = (chunks.q > 0).sum()/len(chunks)
    lim = max(np.abs(chunks.q.min()), np.abs(chunks.q.max()))

    bardata = pd.DataFrame()
    for t in np.arange(0, lim, lim/5):
        total = (np.abs(qs.q)>t).sum()
        bardata = bardata.append({
            't':t,
            'pos':(qs.q>t).sum()/total,
            'neg':(qs.q<-t).sum()/total,
            'total':total}, ignore_index=True)
    ax.bar(bardata.t-lim/15, bardata.pos, lim/15, color='red', linewidth=0.2, label='+')
    ax.bar(bardata.t, bardata.neg, lim/15, color='blue', linewidth=0.2, label='-')
    ax.set_xlim(-lim/15, lim-lim/5+lim/15)
    ax.set_xticks(bardata.t)
    ax.set_xticklabels([
        '{:.1f}'.format(bardata.t.values[0]*1e4),
        '',
        '{:.1f}'.format(bardata.t.values[2]*1e4),
        '',
        '{:.1f}'.format(bardata.t.values[4]*1e4)
        ])
    ax.set_ylim(0, 1)
    ax.set_yticks(np.arange(0,1.1,0.2))
    ax.set_xlabel(r'$t$ $(\times 10^4)$', fontsize=params.labelfontsize)
    ax.set_ylabel(r'% w/strength $>t$', fontsize=params.labelfontsize)
    ax.legend(loc='upper left', fontsize=5, handlelength=0.7, handleheight=0.5, borderpad=0.1,
            labelspacing=0.1, columnspacing=0.1)
    ax.tick_params(**params.tickprops)

    # show manhattan plot
    if c is not None:
        ax = plt.subplot(gs[i,2])
        print('reading sumstats')
        ahat = pd.read_csv('/groups/price/yakir/data/sumstats.hm3/processed/'+pheno+'.sumstats.gz',
                sep='\t').rename(columns={'Z':'Ztrait'})
        print('reading snps')
        snps = pd.read_csv('/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/'+
                '1000G.EUR.QC.'+str(c)+'.bim', header=None,
                names=['CHR','SNP','CM','BP','A1','A2'], delim_whitespace=True)
        snps = snps[(snps.CHR==c)&(snps.BP>=start*1e6)&(snps.BP<=end*1e6)]
        print(snps.BP.min(), snps.BP.max())
        print('merging')
        snps = pd.merge(snps, ahat, on='SNP', how='inner')
        typed = snps.A1_y.notnull()
        snps.loc[typed, 'A1_x'] = snps[typed].A1_y
        snps.loc[typed, 'A2_x'] = snps[typed].A2_y
        snps.rename(columns={'A1_x':'A1', 'A2_x':'A2'}, inplace=True)

        if twas is not None:
            from rpy2.robjects import r, pandas2ri
            import gprim.annotation as ga
            r['load'](twas)
            wgt = r['wgt.matrix']
            top1index = np.where(
                    np.array(r['colnames'](wgt)) == 'top1')[0][0]

            # create dataframe
            mysnps = pd.DataFrame(
                    np.array(r['snps']).T,
                    columns=['CHR','SNP','CM','BP','A1','A2'])
            mysnps['Z'] = np.array(wgt)[:,top1index]
            print(mysnps.BP.min(), mysnps.BP.max())

            #merge
            snps = ga.reconciled_to(snps,
                    mysnps[['SNP','A1','A2','Z']], 'Z')
            ax.scatter(snps.Z, snps.Ztrait, s=1, linewidth=0, c='purple', alpha=0.8)
            ax.set_xlim(1.2*min(snps.Z), 1.2*max(snps.Z))
            ax.set_ylim(1.2*min(snps.Ztrait), 1.2*max(snps.Ztrait))
            ax.set_xticks([z for z in np.arange(-18,18,6) if z > 1.2*min(snps.Z) and z < 1.2*max(snps.Z)])
            ax.set_yticks([z for z in np.arange(-50,50,10) if z > 1.2*min(snps.Ztrait) and z < 1.2*max(snps.Ztrait)])
            ax.set_xlabel(r'$Z_{'+tf+'}$', fontsize=params.labelfontsize)
            ax.set_ylabel(r'$Z_{'+prettypheno+'}$', fontsize=params.labelfontsize)
        else:
            import scipy.stats as st
            ax.locator_params(nbins=2)
            snps['logp'] = -np.log10(st.chi2.sf(snps.Ztrait**2,1))
            sig = snps[snps.logp >= -np.log10(5e-8)]
            nonsig = snps[snps.logp < -np.log10(5e-8)]
            ax.scatter(sig.BP/1e6, sig.logp,
                    s=3, linewidth=0, c='purple')
            ax.scatter(nonsig.BP/1e6, nonsig.logp,
                    s=1, linewidth=0, c='gray', alpha=0.5)
            ax.axhline(y=-np.log10(5e-8), **params.sig_thresh_line_props)
            ax.set_ylim(0, 1.2*max(sig.logp))
            ax.set_xlim(start, end)
            ax.set_xlabel(r'Position (Mb)', fontsize=params.labelfontsize)
            ax.set_ylabel(r'$-\log_{10}(p)$ ('+prettypheno+')',
                    fontsize=params.labelfontsize)
        ax.tick_params(**params.tickprops)

# finish up
sns.despine()
gs.tight_layout(fig)

print('saving figure')
fs.makedir_for_file(outname)
plt.savefig(outname)
plt.close()
