import numpy as np, pandas as pd
from processing_helpers import *
from analysis_helpers import *
from enrichments import *
from mri_maps import *
from gradientVersion import *

def get_gandal_genes(which='rnaseq'):
    disorders = ['ASD', 'MDD', 'SCZ', 'BD']
    if which == 'microarray':

        gandal_genes = (pd.read_csv("../data/gandal_genes_microarray.csv")
                        .rename({'external_gene_id':'gene'},axis=1)
                        .rename({f'{d}.beta_log2FC':f'{d}.log2FC' for d in disorders}, axis=1) 
                        .set_index('gene')
                        .assign(rank = lambda x: x.groupby('gene')['start_position'].rank(method='first'))
                        .loc[lambda x: x['rank']==1]
                        # .loc[lambda x: x.groupby('gene')['rank'].max()>1]
                        # .sort_index()
                       )
    elif which == 'rnaseq':
        gandal_genes = (pd.read_csv("../data/gandal_genes_rnaseq.csv")
                        .rename({'gene_name':'gene'},axis=1) 
                        .rename({f'{d}.fdr':f'{d}.FDR' for d in disorders}, axis=1)                         
                        .set_index('gene')
                        .assign(rank = lambda x: x.groupby('gene')['transcript_length'].rank(method='first'))
                        .loc[lambda x: x['rank']==1]
                        # .loc[lambda x: x.groupby('gene')['rank'].max()>1]
                        # .sort_index()
                       )

    return gandal_genes


def get_gene_corr(weights, null_weights, gandal_genes, sig_thresh=1, disorders = ['ASD', 'MDD', 'SCZ']):
    """
    Correlate gene weights with log2FC
    """
    true_corrs = {}
    null_corrs = {}

    for d in disorders:
        # Skip if disorder not present
        if f'{d}.log2FC' not in gandal_genes.columns:
            continue

        # Don't overwrite the input data in the loop!
        genes_log2FC = gandal_genes.copy()            

        # Filter for significant genes only
        genes_log2FC = genes_log2FC.loc[lambda x: x[f'{d}.FDR'] < sig_thresh, :]            

        # Find matching genes that have non-null differential expression for this disorder
        genes_log2FC = genes_log2FC.loc[lambda x: ~np.isnan(x[f'{d}.log2FC'])]
        genes = set(weights.index).intersection(genes_log2FC.index)
        gene_mask = np.isin(weights.index, list(genes))
        genes_log2FC = genes_log2FC.loc[genes,:]

        # Extract differential expression
        genes_log2FC = genes_log2FC.loc[:, f'{d}.log2FC']

        true = np.corrcoef(x=genes_log2FC.values, y=weights.loc[genes, :].values, rowvar=False)[:1,1:]

        null_corr_disorder = np.zeros((null_weights.shape[2], 3))
        for i in range(3):
            nulls = null_weights[gene_mask, i, :]
            corr = np.corrcoef(x=genes_log2FC.values, y=nulls, rowvar=False)[1:,:1]
            null_corr_disorder[:, i] = corr.squeeze()

        true_corrs[d] = pd.DataFrame(true)
        null_corrs[d] = pd.DataFrame(null_corr_disorder)

    true_corrs = pd.concat(true_corrs).set_axis(['G1','G2','G3'], axis=1).droplevel(1)
    null_corrs = pd.concat(null_corrs).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)

    null_p = compute_null_p(true_corrs, null_corrs, adjust='fdr_bh')
    return null_p


def get_gene_map_corr(version, null_scores, gandal_genes, posneg, method='mean', sig_thresh=1, 
                      disorders = ['ASD','MDD','SCZ'], return_maps=False):
    
    disorders = [d for d in disorders if f'{d}.log2FC' in gandal_genes.columns]

    # Set non sig changes to 0
    for d in disorders:
        gandal_genes[f'{d}.log2FC'] = np.where(gandal_genes[f'{d}.FDR'] >= sig_thresh, 0, gandal_genes[f'{d}.log2FC'])

    gene_log2FC = (gandal_genes
                     .loc[:, [f'{d}.log2FC' for d in disorders]]
                     .rename({f'{d}.log2FC':f'{d}' for d in disorders}, axis=1)
                     .melt(ignore_index=False, var_name='label', value_name='weight')
                     .reset_index()
                    )

    if posneg == 'pos':
        gene_log2FC['weight'] = gene_log2FC['weight'].clip(lower=0)
    elif posneg == 'neg':
        gene_log2FC['weight'] = gene_log2FC['weight'].clip(upper=0)*-1
    elif posneg == 'abs':
        gene_log2FC['weight'] = gene_log2FC['weight'].abs()

    gene_maps = make_gene_maps(version, gene_log2FC, normalize='std', method=method)
    
    if return_maps:
        return gene_maps
    
    true_corrs = get_corrs(version.clean_scores(), gene_maps)
    null_corrs = corr_nulls_from_grads(null_scores, version.clean_scores(), gene_maps)
    null_corrs_p = get_null_p(true_corrs, null_corrs, adjust='fdr_bh').rename({'map':'label'},axis=1)
    
    return null_corrs_p


def get_gene_sig(weights, null_weights, gandal_genes, posneg=None, posneg_weights=None, disorders = ['ASD', 'MDD', 'SCZ']):
    genes_sig_dict = {}
    
    for d in disorders:
        # Skip if disorder not present
        if f'{d}.log2FC' not in gandal_genes.columns:
            continue
        
        # Don't overwrite the input data in the loop!
        gene_sig = gandal_genes.copy()
        
        # Filter for up or down regulated
        if posneg == 'pos':
            gene_sig = gene_sig.loc[lambda x: x[f'{d}.log2FC'] > 0]
        elif posneg == 'neg':
            gene_sig = gene_sig.loc[lambda x: x[f'{d}.log2FC'] < 0]
        elif posneg == 'abs':
            pass
        
        # Filter for significant only
        genes_sig = gene_sig.loc[lambda x: x[f'{d}.FDR'] < .05]
        genes_sig_dict[d] = pd.Series(genes_sig.index)

    # Combine across disorders
    genes_sig = pd.concat(genes_sig_dict).reset_index(0).rename({'gene_name':'gene','level_0':'label'},axis=1)

    true_mean, null_mean = compute_enrichments(weights, null_weights, genes_sig, posneg=posneg_weights)
    null_p = compute_null_p(true_mean, null_mean, adjust='fdr_bh')

    return null_p


def get_gene_dot(weights, null_weights, gandal_genes, posneg=None, posneg_weights=None, sig_thresh=1, disorders = ['ASD', 'MDD', 'SCZ']):
    """
    x
    """
    weights_genes = weights.index
    weights = weights.copy().values
    nulls = null_weights.copy()
    
    if posneg_weights =='abs':
        weights = np.abs(weights)
        nulls = np.abs(nulls)
    elif posneg_weights=='pos':
        weights = np.where(weights<0, 0, weights)
        nulls = np.where(nulls<0, 0, nulls)
    elif posneg_weights=='neg':
        weights = np.where(weights>0, 0, weights)
        nulls = np.where(nulls>0, 0, nulls)
    
    true_dot = {}
    null_dot = {}
    for d in disorders:
        # Skip if disorder not present
        if f'{d}.log2FC' not in gandal_genes.columns:
            continue

        # Don't overwrite the input data in the loop!
        genes_log2FC = gandal_genes.copy()           
            
        # Filter for significant genes only
        genes_log2FC = genes_log2FC.loc[lambda x: x[f'{d}.FDR'] < sig_thresh, :]            
            
        # Find matching genes that have non-null differential expression for this disorder
        genes_log2FC = genes_log2FC.loc[lambda x: ~np.isnan(x[f'{d}.log2FC'])]
        matched_genes = set(weights_genes).intersection(genes_log2FC.index)
        gene_mask = np.isin(weights_genes, list(matched_genes))
        
        # Extract differential expression
        genes_log2FC = genes_log2FC.loc[matched_genes, f'{d}.log2FC']
        if posneg == 'pos':
            genes_log2FC = genes_log2FC.clip(lower=0)
        elif posneg == 'neg':
            genes_log2FC = genes_log2FC.clip(upper=0)*-1
        elif posneg == 'abs':
            genes_log2FC = genes_log2FC.abs()

        # Filter weights for matching genes
        weights_matched = weights[gene_mask, :3]
        null_weights_matched = np.swapaxes(nulls[gene_mask, :, :], 0, 2)
                
        # Compute dot product
        true_dot[d] = pd.Series(weights_matched.T @ genes_log2FC)
        null_dot[d] = pd.DataFrame(np.dot(null_weights_matched, genes_log2FC))

    true_dot = pd.concat(true_dot).unstack(1).set_axis(['G1','G2','G3'], axis=1)
    null_dot = pd.concat(null_dot).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)
    
    null_p = compute_null_p(true_dot, null_dot, adjust='fdr_bh')
    
    return null_p


