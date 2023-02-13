import numpy as np, pandas as pd
from processing_helpers import *
from analysis_helpers import *
from enrichments import *
from mri_maps import *
from gradientVersion import *

def get_disgenet_genes():
    disgenet_genes = (pd.read_csv("../data/disgenet_genes.csv")
        .assign(score_pct = lambda x: x.groupby('Disease')['Score_gda'].rank(pct=True))
        # .loc[lambda x: x['score_pct'] >= 0.75]
        .rename({'Disease':'label','Gene':'gene'},axis=1)
        .replace({
            'Autism Spectrum Disorders':'ASD',
            'Major Depressive Disorder':'MDD',
            'Schizophrenia':'SCZ'
        })
    )
    return disgenet_genes 

def get_gandal_genes(which='microarray', disorders = ['ASD', 'SCZ', 'MDD'], sig_level=.05):

    replace_dict = {'01-Mar':'MARCH1', '02-Mar':'MARCH2', '03-Mar':'MARCH3', '04-Mar':'MARCH4', '05-Mar':'MARCH5', '06-Mar':'MARCH6', '07-Mar':'MARCH7', '08-Mar':'MARCH8', 
                    '09-Mar':'MARCH9', '10-Mar':'MARCH10', '11-Mar':'MARCH11',
                    '01-Sep':'SEPT1', '02-Sep':'SEPT2', '03-Sep':'SEPT3', '04-Sep':'SEPT4', '05-Sep':'SEPT5', '06-Sep':'SEPT6', '07-Sep':'SEPT7', '08-Sep':'SEPT8',
                    '09-Sep':'SEPT9', '10-Sep':'SEPT10', '11-Sep':'SEPT11', '12-Sep':'SEPT12', '13-Sep':'SEPT13', '14-Sep':'SEPT14', '15-Sep':'SEPT15', 
                    '01-Dec':'DECR1', '02-Dec':'DECR2'}

    if which == 'microarray':
        gandal_genes = (pd.read_csv("../data/gandal_genes_microarray.csv")
                        .rename({'external_gene_id':'gene'},axis=1)
                        .replace({'gene': replace_dict})
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
                        .replace({'gene':replace_dict})
                        .rename({f'{d}.fdr':f'{d}.FDR' for d in disorders}, axis=1)                         
                        .set_index('gene')
                        .assign(rank = lambda x: x.groupby('gene')['transcript_length'].rank(method='first'))
                        .loc[lambda x: x['rank']==1]
                        # .loc[lambda x: x.groupby('gene')['rank'].max()>1]
                        # .sort_index()
                       )
        
    for d in disorders:
        gandal_genes[f'{d}.sig'] = gandal_genes[f'{d}.FDR'] < sig_level

    gandal_sig = (gandal_genes
            .loc[:,[f'{d}.sig' for d in disorders]]
            .melt(ignore_index=False, var_name='label', value_name='sig')
            .assign(label = lambda x: x['label'].str.replace('.sig','', regex=False))
            .set_index('label', append=True)
            )

    gandal_log2FC = (gandal_genes
            .loc[:,[f'{d}.log2FC' for d in disorders]]
            .melt(ignore_index=False, var_name='label', value_name='log2FC')
            .assign(label = lambda x: x['label'].str.replace('.log2FC','', regex=False))
            .set_index('label', append=True)
            .join(gandal_sig)
        )

    return gandal_log2FC


def get_gene_weights_with_labels(weights, gandal_genes, 
                                 disorders = ['ASD','SCZ','MDD']):
    gandal_fdr = (gandal_genes
                  .loc[:,[f'{d}.FDR' for d in disorders]]
                  .melt(ignore_index=False, var_name='disorder', value_name='FDR')
                  # .loc[lambda x: x['FDR']]
                  .assign(disorder = lambda x: 
                          x['disorder'].str.replace('.FDR','', regex=False))
                  .set_index('disorder', append=True)
                 )

    gandal_log2FC = (gandal_genes
                  .loc[:,[f'{d}.log2FC' for d in disorders]]
                  .melt(ignore_index=False, var_name='disorder', value_name='log2FC')
                  .assign(disorder = lambda x: 
                          x['disorder'].str.replace('.log2FC','', regex=False))
                  .set_index('disorder', append=True)
                  .join(gandal_fdr)
                  .reset_index('disorder')
                  .loc[lambda x: x['FDR'] < 0.05]
                 )

    weights_labels = (weights
                      .join(gandal_log2FC)     
                      .rename({'disorder':'label'}, axis=1)                   
                      .fillna({'label':'none'})
                      .assign(cluster = lambda x: 
                              np.select([x['log2FC']>0,x['log2FC']<0],
                                        ['upregulated','downregulated'], 
                                        default=None)
                             )
                     )
    return weights_labels


def get_rare_genes():
    rare_genes = pd.read_csv("../data/rare_genes_konrad.csv").melt(var_name='label',value_name='gene').dropna()
    return rare_genes


def np_pearson_corr(x, y):
    """
    Fast pearson correlation (don't cross-correlate all the nulls with each other)
    """
    xv = x - x.mean(axis=0)
    yv = y - y.mean(axis=0)
    xvss = (xv * xv).sum(axis=0)
    yvss = (yv * yv).sum(axis=0)
    result = np.matmul(xv.transpose(), yv) / np.sqrt(np.outer(xvss, yvss))
    # bound the values to -1 to 1 in the event of precision issues
    return np.maximum(np.minimum(result, 1.0), -1.0)

def get_gene_corr(weights, null_weights, gandal_genes, sig_thresh=1, adjust='fdr_bh', disorders = ['ASD', 'MDD', 'SCZ']):
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

        # Take only significant genes
        genes_log2FC = genes_log2FC.loc[lambda x: x[f'{d}.FDR']<sig_thresh, f'{d}.log2FC']
        
        # Find matching genes and filter differential expression
        genes_matched = set(weights.index).intersection(genes_log2FC.index)
        genes_log2FC = genes_log2FC.loc[genes_matched].values
        # Also filter true and null weights
        # np.searchsorted to find indices of matched genes in same order as genes_matched
        genes_idxs = np.searchsorted(weights.index, list(genes_matched))
        true = weights.values[genes_idxs, :]
        nulls = null_weights[genes_idxs, :, :]
        
        # For each PC, get the correlation of the true weights and all the nulls
        true_corr_disorder = np.zeros(3)
        null_corr_disorder = np.zeros((null_weights.shape[2], 3))
        for i in range(3):
            true_corr_disorder[i] = np_pearson_corr(genes_log2FC, true[:,i]).squeeze()
            null_corr_disorder[:, i] = np_pearson_corr(genes_log2FC, nulls[:,i,:]).squeeze()
        
        true_corrs[d] = pd.DataFrame(true_corr_disorder).T
        null_corrs[d] = pd.DataFrame(null_corr_disorder) 

    true_corrs = pd.concat(true_corrs).set_axis(['G1','G2','G3'], axis=1).droplevel(1)
    null_corrs = pd.concat(null_corrs).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)

    null_p = compute_null_p(true_corrs, null_corrs, adjust=adjust)
    return null_p



def get_gene_corr2(weights, null_weights, genes_log2FC, group='disorder', sig_thresh=1, adjust='fdr_bh'):
    """
    Correlate gene weights with log2FC
    """
    true_corrs = {}
    null_corrs = {}

    for g in genes_log2FC[group].unique():
        # Don't overwrite the input data in the loop!
        _genes_log2FC = genes_log2FC.copy().loc[lambda x: x[group]==g]

        # Take only significant genes
        _genes_log2FC = _genes_log2FC.loc[lambda x: x['FDR']<sig_thresh, 'log2FC']
        
        # Merge duplicates
        _genes_log2FC = _genes_log2FC.groupby('gene').mean()

        # Find matching genes and filter differential expression
        genes_matched = list(set(weights.index).intersection(_genes_log2FC.index))
        _genes_log2FC = _genes_log2FC.loc[genes_matched].values
        # Also filter true and null weights
        # np.searchsorted to find indices of matched genes in same order as genes_matched
        genes_idxs = np.searchsorted(weights.index, list(genes_matched))
        true = weights.values[genes_idxs, :]
        nulls = null_weights[genes_idxs, :, :]
        
        # For each PC, get the correlation of the true weights and all the nulls
        true_corr_group = np.zeros(3)
        null_corr_group = np.zeros((null_weights.shape[2], 3))
        for i in range(3):
            true_corr_group[i] = np_pearson_corr(_genes_log2FC, true[:,i]).squeeze()
            null_corr_group[:, i] = np_pearson_corr(_genes_log2FC, nulls[:,i,:]).squeeze()
        
        true_corrs[g] = pd.DataFrame(true_corr_group).T
        null_corrs[g] = pd.DataFrame(null_corr_group) 

    true_corrs = pd.concat(true_corrs).set_axis(['G1','G2','G3'], axis=1).droplevel(1)
    null_corrs = pd.concat(null_corrs).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)

    # return true_corrs, null_corrs
    null_p = compute_null_p(true_corrs, null_corrs, adjust=adjust)
    return null_p




def get_gene_map_corr(version, null_scores, gandal_genes, posneg, method='mean', sig_thresh=1, 
                      disorders = ['ASD','MDD','SCZ'], return_maps=False):
    
    disorders = [d for d in disorders if f'{d}.log2FC' in gandal_genes.columns]

    gene_log2FC = gandal_genes.copy()
    
    # Set non sig changes to 0
    for d in disorders:
        gene_log2FC[f'{d}.log2FC'] = np.where(gene_log2FC[f'{d}.FDR'] >= sig_thresh, 0, gene_log2FC[f'{d}.log2FC'])

    gene_log2FC = (gene_log2FC
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


def get_gene_sig(weights, null_weights, gandal_genes, posneg=None, posneg_weights=None, disorders = ['ASD', 'MDD', 'SCZ'], combine=None):
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
    
    # Take superset if desired
    if combine=='union':
        genes_sig['label'] = '-'.join(genes_sig['label'].unique())
        genes_sig = genes_sig.drop_duplicates()
    if combine=='ASD-SCZ intersection':
        ASD_genes = genes_sig.loc[lambda x: x['label'] == 'ASD', 'gene']
        SCZ_genes = genes_sig.loc[lambda x: x['label'] == 'SCZ', 'gene']
        intersection = set(ASD_genes).intersection(SCZ_genes)
        genes_sig = genes_sig.set_index('gene').loc[intersection, :].assign(label = 'ASD-SCZ\nintersection').reset_index()
    
    true_mean, null_mean = compute_enrichments(weights, null_weights, genes_sig, posneg=posneg_weights)
    null_p = compute_null_p(true_mean, null_mean, adjust='fdr_bh')

    return null_p


def get_gene_dot(weights, null_weights, gandal_genes, posneg=None, posneg_weights=None, sig_thresh=1, disorders = ['ASD', 'MDD', 'SCZ']):
    """
    x
    """
    weights = weights.copy()
    null_weights = null_weights.copy()
    if posneg_weights =='abs':
        weights = np.abs(weights)
        null_weights = np.abs(null_weights)
    elif posneg_weights=='pos':
        weights = pd.DataFrame(np.where(weights<0, 0, weights), index=weights.index)
        null_weights = np.where(null_weights<0, 0, null_weights)
    elif posneg_weights=='neg':
        weights = pd.DataFrame(np.where(weights>0, 0, weights), index=weights.index)
        null_weights = np.where(null_weights>0, 0, null_weights)    
    
    true_dot = {}
    null_dot = {}
    for d in disorders:
        # Skip if disorder not present
        if f'{d}.log2FC' not in gandal_genes.columns:
            continue          
            
        # Don't overwrite the input data in the loop!
        genes_log2FC = gandal_genes.copy()           
            
        # Take only significant genes
        genes_log2FC = genes_log2FC.loc[lambda x: x[f'{d}.FDR']<sig_thresh, f'{d}.log2FC']

        # Find matching genes and filter differential expression
        genes_matched = set(weights.index).intersection(genes_log2FC.index)
        genes_log2FC = genes_log2FC.loc[genes_matched].values
        # Also filter true and null weights
        # np.searchsorted to find indices of matched genes in same order as genes_matched
        genes_idxs = np.searchsorted(weights.index, list(genes_matched))
        true = weights.values[genes_idxs, :]
        nulls = null_weights[genes_idxs, :, :]
        
        # Clip or take absolute
        if posneg == 'pos':
            genes_log2FC = genes_log2FC.clip(min=0)
        elif posneg == 'neg':
            genes_log2FC = genes_log2FC.clip(max=0)*-1
        elif posneg == 'abs':
            genes_log2FC = genes_log2FC.abs()

        # Swap nulls axes for dot product
        nulls = np.swapaxes(nulls[:, :, :], 0, 2)
                
        # Compute dot product
        true_dot[d] = pd.Series(true.T @ genes_log2FC)
        null_dot[d] = pd.DataFrame(np.dot(nulls, genes_log2FC))

    true_dot = pd.concat(true_dot).unstack(1).set_axis(['G1','G2','G3'], axis=1)
    null_dot = pd.concat(null_dot).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)
    
    null_p = compute_null_p(true_dot, null_dot, adjust='fdr_bh')
    
    return null_p


