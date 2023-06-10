# Functions for enrichment analysis

import numpy as np, pandas as pd
from scipy.stats import percentileofscore
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import gseapy as gp
from processing_helpers import *


def shuffle_gene_weights(weights, n=100, rank=False):#, n_components=3):
    """
    Make null model by randomizing gene weights / ranks
    """
    n_components=weights.shape[1]
    null_weights = np.repeat(weights.values[:,:n_components, np.newaxis], n, axis=2)
    # null_weights = np.take_along_axis(null_weights, np.random.randn(*null_weights.shape).argsort(axis=0), axis=0)
    for c in range(n_components):
        for i in range(n):
            np.random.shuffle(null_weights[:,c,i])
    
    if rank:
        null_ranks = null_weights.argsort(axis=0).argsort(axis=0)
        return null_ranks
    else:
        return null_weights


def match_genes(gene_labels, weights):
    """
    Make mask of which genes from each gene list are in the PCs
    """
    genes = weights.index
    gene_masks = {}
    gene_counts = {}
    for l in gene_labels['label'].unique():
        genes_in_label = gene_labels.query("label == @l")['gene']
        matches = np.isin(genes, genes_in_label)
        if sum(matches)>0:
            gene_masks[l] = matches
        
        gene_counts[l] = pd.Series({
            'n_genes': len(genes_in_label),
            'n_matches': sum(matches)
        })
    gene_counts = pd.concat(gene_counts).unstack()
    
    return gene_masks, gene_counts


def compute_enrichments(weights, null_weights, gene_labels, 
                        how='mean', norm=False, posneg=None):
    """
    Compute scores for each gene label, either mean, or median rank
    """
    n_components = weights.shape[1]
    axis_names = list(weights.columns)
    gene_masks, gene_counts = match_genes(gene_labels, weights)
    
    weights = weights.copy().values
    nulls = null_weights.copy()
    # Take absolute values of standardized weights
    if norm:
        weights = StandardScaler().fit_transform(weights)
        for i in range(nulls.shape[2]):
            nulls[:,:,i] = StandardScaler().fit_transform(nulls[:,:,i])
            
    if posneg =='abs':
        weights = np.abs(weights)
        nulls = np.abs(nulls)
    elif posneg=='pos':
        weights = np.where(weights<0, np.nan, weights)
        nulls = np.where(nulls<0, np.nan, nulls)
    elif posneg=='neg':
        weights = np.where(weights>0, np.nan, weights)
        nulls = np.where(nulls>0, np.nan, nulls)

    true_enrichments = {}
    null_enrichments = {}
    
    for label, mask in gene_masks.items():
        if how == 'mean':
            true_enrichments[label] = pd.Series(np.nanmean(weights[mask, :], axis=0))
            null_enrichments[label] = pd.DataFrame(np.nanmean(nulls[mask, :, :], axis=0)).T
        elif how == 'median': #### not working
            true_ranks = weights.argsort(0).argsort(0)
            true_enrichments[label] = pd.Series(np.nanmedian(true_ranks[mask, :], axis=0))
            null_enrichments[label] = pd.DataFrame(np.nanmedian(nulls[mask, :, :], axis=0)).T

    true_enrichments = pd.concat(true_enrichments).unstack(1).set_axis(axis_names, axis=1)
    null_enrichments = pd.concat(null_enrichments).set_axis(axis_names, axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)

    return true_enrichments, null_enrichments, gene_counts
    
    
def compute_null_p(true_enrichments, null_enrichments, 
                   gene_counts=None, adjust='fdr_bh', adjust_by_label=False,
                   order=None):
    """
    Compute null p values
    """
    null_pct = np.zeros(true_enrichments.shape)
    for m, label in enumerate(true_enrichments.index):
        for i in range(true_enrichments.shape[1]):
            nulls_ = null_enrichments.set_index('label').loc[label].iloc[:,i]
            true_ = true_enrichments.iloc[m, i]
            pct = percentileofscore(nulls_, true_)/100
            null_pct[m, i] = pct
            
    true_mean = true_enrichments.stack().rename('true_mean')

    null_mean = (null_enrichments
                 .groupby('label').agg(['mean','std'])
                 .stack(0)
                 .rename_axis([None,None])
                 .set_axis(['null_mean', 'null_std'], axis=1)
                )
            
    null_p = (pd.DataFrame(null_pct,
                           index=true_enrichments.index,
                           columns=true_enrichments.columns)
              .stack().rename('pct').to_frame()
              .join(true_mean)
              .join(null_mean)
              .assign(z = lambda x: (x['true_mean'] - x['null_mean'])/x['null_std'])
              .assign(pos = lambda x: [pct > 0.5 for pct in x['pct']])
              .assign(p = lambda x: [(1-pct)*2 if pct>0.5 else pct*2 for pct in x['pct']]) # x2 for bidirectional
             )
    
    # Apply multiple comparisons
    if adjust is not None:
        # Adjust across axes only (not by label)?
        if adjust_by_label:
            null_p = (null_p
                .assign(q = lambda x: x.groupby(level=0)
                                       .apply(lambda y: pd.Series(multipletests(y['p'], method=adjust)[1], index=y.index))
                                       .reset_index(0, drop=True) # make index match
                                       )
                .assign(sig = lambda x: x['q'] < .05)
                # .assign(q_abs = lambda x: [1-q if pos else q for pos, q in zip(x['pos'], x['q'])])
                )
        else:
            null_p = (null_p
                .assign(q = lambda x: multipletests(x['p'], method=adjust)[1])
                .assign(sig = lambda x: x['q'] < .05)
                # .assign(q_abs = lambda x: [1-q if pos else q for pos, q in zip(x['pos'], x['q'])])
                )
    else:
        null_p = (null_p
             .assign(q = lambda x: x['p'])
             .assign(sig = lambda x: x['q'] < .05)
            )
    
    null_p = (null_p
              .reset_index()
              .rename({'level_0':'label', 'level_1':'C'}, axis=1)
             )
    
    if gene_counts is not None:
        null_p = null_p.join(gene_counts, on='label')
    
    # Fix order of gene labels
    if order is None:
        order = true_enrichments.index
    null_p = (null_p
              .assign(label = lambda x: pd.Categorical(x['label'], ordered=True, categories=order))
              .sort_values('label')
         )

    return null_p



def make_axis_quantiles(weights, q=10, labels=None, levels=None, na_value=np.NaN):
    """
    Split weights by axis quantiles, then match to labels
    """
    weights_with_quantiles = (weights
            .melt(ignore_index=False, var_name='C', value_name='C_score')
            .reset_index().rename({'index':'gene'}, axis=1)
            .assign(C_quantile = lambda x: x.groupby('C').apply(lambda y: pd.qcut(y['C_score'], q=q, labels=range(q))).reset_index(0, drop=True))
    )
    if labels is None:
        return weights_with_quantiles
    else:
        if levels is None:
            levels = labels['label'].unique()
        return (weights_with_quantiles
                .join(labels.set_index('gene'), on='gene')
                .fillna({'label': na_value})
                .assign(label = lambda x: pd.Categorical(x['label'], ordered=True, categories=levels))
                .sort_values(['label','C_score'], ascending=False)
                .assign(rank_in_quantile = lambda x: x.groupby(['C','C_quantile']).cumcount()+1)
                .sort_values(['C_quantile','rank_in_quantile'])
            )


def make_label_quantiles_by_axis(weights, labels, q=10):
    """
    Split labels into quantiles by axis weights
    """
    label_percentiles = (labels.loc[:, ['label','gene']]
        .join(weights.melt(ignore_index=False, var_name='C',value_name='C_score'), on='gene')
        .dropna()
        .reset_index(drop=True)
        .assign(C_quantile = lambda x: x.groupby(['C','label'])['C_score'].apply(lambda y: pd.qcut(y, q=q, labels=range(1,q+1
        ))))
    )
    return label_percentiles



def make_gene_maps(version, gene_labels, atlas='hcp', normalize='std', method='mean', order=None, use_weights=True):
    gene_maps = {}
    for label in gene_labels['label'].unique():
        gene_labels_ = gene_labels.set_index('label').loc[label]
        if len(gene_labels_)== 1: 
            gene_labels_ = gene_labels_.to_frame().T
        mask = np.intersect1d(gene_labels_['gene'], version.expression.columns)
        if len(mask) == 0:
            continue
        gene_label_expression_ = version.expression.loc[:, mask]

        if 'weight' in gene_labels.columns and use_weights:
            gene_label_weights = gene_labels_.set_index('gene').loc[mask, 'weight']
            gene_label_expression_ = gene_label_expression_ * gene_label_weights 
        
        # Mean expression
        gene_mean = gene_label_expression_.mean(axis=1)
        gene_eigengene = PCA(n_components=1).fit_transform(gene_label_expression_.dropna(axis=1)).squeeze()
        gene_eigengene = pd.Series(gene_eigengene, index=gene_label_expression_.index)
        corr = gene_mean.corr(gene_eigengene)
        # print(f"label={label}: corr(mean expression, eigengene) = {corr}")
        if corr < 0:
            gene_eigengene *= -1
        
        if method=='mean':
            gene_maps[label] = gene_mean
        elif method=='eigengene':
            gene_maps[label] = gene_eigengene

    gene_maps = pd.concat(gene_maps, axis=1)
    
    if order is not None:
        gene_maps = gene_maps.loc[:,order]
    
    if normalize == 'std':
        gene_maps = gene_maps.apply(lambda x: (x-np.mean(x))/np.std(x))
    elif normalize == 'fisher':
        gene_maps = gene_maps.apply(lambda x: 0.5*(1 + (np.exp(x) - np.exp(-x)) / (np.exp(x) + np.exp(-x)) ))
    
    if atlas == 'hcp':
        gene_maps = gene_maps.join(get_labels_hcp())#.set_index('label')
    elif atlas == 'dk':
        gene_maps = gene_maps.join(get_labels_dk())#.set_index('label')
        
    return gene_maps



def make_gene_corrs(genes, version, maps):
    genes = set(genes).intersection(version.expression.columns)
    gene_maps = (version.expression.loc[:,genes]
                 .join(get_labels_hcp()).set_index('label'))

    corrs = (pd.concat([gene_maps, maps], axis=1)
             .corr().iloc[:len(genes),len(genes):]
             .reset_index(drop=True)
            )
    
    for col in corrs:
        corrs[col] = corrs[col].sort_values(ignore_index=True)
        
    return corrs


### GO enrichment

def get_target_enrichments(target, background,
                                      p_threshold=0.05,
                                      FDR_threshold=1
                                     ):
    results = gp.enrichr(gene_list=target,
           background=background,
           gene_sets="../data/c5.go.bp.v7.5.1.symbols.gmt",
           organism='human',
           cutoff=0.05,
           outdir='../outputs/enrichr',
           no_plot=True,
           verbose=True).results
    
    results = (results
               .loc[:, ['Term', 'Overlap', 'P-value', 'Adjusted P-value']]
               .set_axis(['term', 'overlap', 'p', 'fdr'], axis=1)
               .sort_values('p')
               .loc[lambda x: x['p'] < p_threshold]
               .loc[lambda x: x['fdr'] < FDR_threshold]
               .assign(term = lambda x: x['term'].str[5:].str.lower().str.replace('_',' '))
              )
    
    return results