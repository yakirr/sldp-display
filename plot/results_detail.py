from __future__ import print_function, division
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
import matplotlib.cm as cm
import scipy.stats as st
# from rpy2.robjects import r, pandas2ri
import statutils.vis as vis
import gprim.annotation as ga
from plot import params

def plot_ahat_vs_Rv(ax, indir, pheno, annot, (xmin, xmax, xexp), (ymin, ymax, yexp)):
    annotname = annot+',diff,'+annot+',-1.R'
    print(pheno, annot)
    print(indir+pheno+'.'+annotname+'.toplot.gz')
    data = pd.read_csv(indir+pheno+'.'+annotname+'.toplot.gz', sep='\t')
    chunks = pd.read_csv(indir+pheno+'.'+annotname+'.chunks', sep='\t')
    mu = chunks.q.sum()/chunks.r.sum()
    xmult = np.power(10, xexp)
    ymult = np.power(10, yexp)

    # make plot of ahat vs Rv
    ax.set_xticks([0])
    ax.set_yticks([0])
    xs, ys = vis.scatter_b(
            data.Rv_resid*xmult,
            data.ahat_resid*ymult,
            binsize=4000, extreme_only=5, errorbars=False,
            ax=ax,
            color='purple',
            fmt='o',
            linewidth=0.5,
            markersize=2)
    print(min(xs), max(xs), min(ys), max(ys))
    ax.plot([min(xs), max(xs)], [mu*min(xs)*ymult/xmult, mu*max(xs)*ymult/xmult],
        c='gray', linewidth=0.8, dashes=[2,2])
    ax.set_xlim(xmin, xmax); ax.set_xticks([xmin, xmax])
    ax.set_ylim(ymin, ymax); ax.set_yticks([ymin, ymax])
    ax.set_xlabel(r'avg $Rv$ $(\times 10^{'+str(-xexp)+r'})$',
            fontsize=params.labelfontsize)
    ax.set_ylabel(r'avg $\hat\alpha$ $(\times 10^{'+str(-yexp)+r'})$',
            fontsize=params.labelfontsize)
    ax.tick_params(**params.tickprops)

    numbers = pd.DataFrame()
    numbers['mean alphahat'] = xs / xmult
    numbers['mean Rv'] = ys / ymult
    return numbers

def manhattan(fig, subplotspec, pheno, tf, c, start, end, gstart, gend,
        ytick=None, yrange=None,
        twas=None, show_gene_loc=False):
    def add_twas(snps, twas):
        # read TWAS weights
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
        return ga.reconciled_to(snps,
                mysnps[['SNP','A1','A2','Z']], 'Z')

    # read trait sumstats and process
    print('reading sumstats')
    ahat = pd.read_csv(params.sumstats+'processed/'+pheno+'.sumstats.gz',
            sep='\t').rename(columns={'Z':'Ztrait'})
    print('reading snps')
    snps = pd.read_csv(params.kg3+'plink_files/'+
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
    print('computing logp etc')
    snps['logp'] = -np.log10(st.chi2.sf(snps.Ztrait**2,1))
    snps['sig'] = (snps.logp >= -np.log10(5e-8))
    snps['pol_logp'] = snps.logp * np.sign(snps.Ztrait)
    print(snps.pol_logp.min(), snps.pol_logp.max(), snps.logp.min(), snps.logp.max())

    if twas is None:
        # get axis
        ax = plt.subplot(subplotspec)

        # plot
        ax.scatter(snps[snps.sig].BP/1e6, snps[snps.sig].logp,
                s=4, linewidth=0, c='purple')
        ax.scatter(snps[~snps.sig].BP/1e6, snps[~snps.sig].logp,
                s=1, linewidth=0, c='gray', alpha=0.5)

        # significance line
        ax.axhline(y=-np.log10(5e-8), **params.sig_thresh_line_props)

        # format axes
        if yrange is None:
            ax.set_ylim(0, 1.2*max(snps.logp))
        else:
            ax.set_ylim(0, yrange)
        if ytick is not None:
            ax.set_yticks([0, ytick])
        ax.set_ylabel(r'$-\log_{10}(p)$',
                fontsize=params.labelfontsize)
    else:
        # make axes for plot and colorbar
        inner_grid = gridspec.GridSpecFromSubplotSpec(1,100, subplotspec)
        ax = plt.subplot(inner_grid[0,:90])
        cax = plt.subplot(inner_grid[0,93:])

        # add twas data and polarize plot
        snps = add_twas(snps, twas)

        # add size information
        snps.loc[snps.sig,'s'] = 4
        snps.loc[~snps.sig,'s'] = 1

        # plot
        plot = snps; # plot = snps[snps.sig]
        cb = ax.scatter(plot.BP/1e6, plot.pol_logp,
                c=np.sign(plot.Z)*np.sqrt(np.abs(plot.Z)), cmap=cm.bwr,
                s=plot.s, linewidth=0)
        # plot = snps[~snps.sig]
        # ax.scatter(plot.BP/1e6, plot.pol_logp,
        #         c='gray', s=plot.s, linewidth=0)

        # colorbar properties
        cb = fig.colorbar(cb, cax=cax, ticks=[]) # could set to [-4.5, 4.5]
        cb.set_label(label=r'$\it{'+tf+'}$ eQTL z-score', size=params.labelfontsize-2)
        cb.outline.set_linewidth(0.05)
        cax.tick_params(**params.tickprops)

        # significance lines and line at y=0
        ax.axhline(y=-np.log10(5e-8), **params.sig_thresh_line_props)
        ax.axhline(y=+np.log10(5e-8), **params.sig_thresh_line_props)
        ax.axhline(y=0, color='black', linewidth=0.5)

        # format axes
        if yrange is None:
            ax.set_ylim(-1.1*max(snps.logp), 1.1*max(snps.logp))
        else:
            ax.set_ylim(-yrange/2, yrange/2)
        if ytick is not None:
            ax.set_yticks([-ytick, 0, ytick])
        ax.set_ylabel(r'polarized $-\log_{10}(p)$',
                fontsize=params.labelfontsize,
                labelpad=-1)
        ax.spines['bottom'].set_color('none')

    # add vertical lines at gene end point if required
    if show_gene_loc:
        ax.axvline(x=gstart, color='green', linewidth=0.5)
        ax.axvline(x=gend, color='green', linewidth=0.5)

    # apply general axis formatting for both types of plots
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.set_xticks([start, end])
    ax.set_xlim(start, end)
    ax.set_xlabel(r'chr {} (Mb)'.format(c), fontsize=params.labelfontsize)
    ax.tick_params(**params.tickprops)

    if twas is None:
        return snps[['CHR','BP','SNP','A1','A2','logp']].rename(
                columns={'logp':'-log10(p)'})
    else:
        return snps[['CHR','BP','SNP','A1','A2','pol_logp','Z']].rename(
                columns={'pol_logp':'polarized -log10(p)',
                    'Z':'IRF1 expr Z-score'})

def plot_q(ax, indir, pheno, annot):
    annotname = annot+',diff,'+annot+',-1.R'
    print(pheno, annot)
    print(indir+pheno+'.'+annotname+'.toplot.gz')
    chunks = pd.read_csv(indir+pheno+'.'+annotname+'.chunks', sep='\t')
    numbers = pd.DataFrame()
    numbers['slope'] = np.sort(chunks.q.values)
    half = np.argmin(numbers.slope**2)
    percent = (numbers.slope > 0).sum()/len(numbers)
    lim = max(np.abs(numbers.slope.min()), np.abs(numbers.slope.max()))

    bardata = pd.DataFrame()
    for t in np.arange(0, lim, lim/5):
        total = (np.abs(numbers.slope)>t).sum()
        bardata = bardata.append({
            't':t,
            'pos':(numbers.slope>t).sum()/total,
            'neg':(numbers.slope<-t).sum()/total,
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
    ax.set_xlabel(r'$t$ $(\times 10^4)$', fontsize=params.labelfontsize-1)
    ax.set_ylabel(r'% w/strength $>t$', fontsize=params.labelfontsize-1)
    ax.legend(loc='upper left', fontsize=5, handlelength=0.7, handleheight=0.5, borderpad=0.1,
            labelspacing=0.1, columnspacing=0.1)
    ax.tick_params(**dict(params.tickprops, labelsize=6))

def enrichment(ax, enrichments, pheno, enrichment_tfs, ticks, num_enrichments=2):
    print('reading enrichments and filtering')
    e = pd.read_csv(enrichments, sep='\t')
    e = e[e.target.str.contains(pheno)]
    containstfs = np.any([e.target.str.contains(tf) for tf in enrichment_tfs], axis=0)
    e = e[containstfs]
    e['x'] = e.mean_in_wgt / e.mean_out_wgt
    e['xabs'] = e.mean_in_wgt - e.mean_out_wgt
    # sortkey=['p_wgt','xabs']
    sortkey=['q','xabs']
    print('sorting by', sortkey)
    top = e.sort_values(by=sortkey)[:2]
    if len(top.geneset.unique()) == 1:
        top = e.sort_values(by=sortkey).iloc[[0,2]]

    print(top[['target','p_wgt','mean_in_wgt','mean_out_wgt']])
    print()
    print([(g,z) for g,z in zip(top.geneset, top.x)])
    print()

    ind = np.arange(2)
    width=0.7
    rects1 = ax.barh(ind,
            top.mean_in_wgt / top.mean_out_wgt, width,
            xerr=top.std_in_wgt / top.mean_out_wgt,
            color='IndianRed',
            error_kw=dict(ecolor='black', lw=1, capsize=0, capthick=0),
            label='in gene set')
    ax.axvline(x=1, **params.enrichment_thresh_line_props)

    ax.margins(y=0.2)
    maxy = max(top.mean_in_wgt / top.mean_out_wgt + top.std_in_wgt / top.mean_out_wgt)
    ax.set_xlim((0, max(ticks)))
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(t)+'x' for t in ticks])
    ax.set_xlabel('SLDP enrichment\nof gene set', fontsize=params.labelfontsize)
    ax.set_yticks([])
    ax.tick_params(**params.tickprops)

    top['Trait'] = top.target.str.split('.').str.get(0)
    top['Annotation'] = top.target.str.split('.').str.get(1)
    top.loc[top.p_wgt == 0, 'p_wgt'] = '<1e-6'
    return top[['Trait','Annotation','geneset','mean_in_wgt','mean_out_wgt',
        'std_in_wgt','std_out_wgt','p_wgt']].rename(
                columns={
                    'geneset':'Gene set',
                    'mean_in_wgt':'Avg covariance in set',
                    'mean_out_wgt':'Avg covariance outside set',
                    'std_in_wgt':'Std error(Avg in set)',
                    'std_out_wgt':'Std error(Avg outside set)',
                    'p_wgt':'p'})
