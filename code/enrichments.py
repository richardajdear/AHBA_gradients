# Functions for enrichment analysis

import numpy as np, pandas as pd
from scipy.stats import percentileofscore

def get_cell_genes(path="../data/jakob_cell_genes.csv"):
    """
    Read cell genes table from Jakob's paper
    """
    cell_genes = pd.read_csv(path)
    
    cell_genes_merged = (cell_genes
     .loc[lambda x: x['Paper'] != 'Lake', :]
     # .drop('Class', axis=1).assign(Class = lambda x: x['Type'])
     .replace({'Class':{'Neuro-Ex':'Neuro', 'Neuro-In':'Neuro'}})
     .melt(id_vars=['Type', 'Paper', 'Cluster', 'Class'], value_name='Gene')
     .loc[lambda x: x['Gene'].notnull(), ['Class', 'Gene']]
     .drop_duplicates()
     # .loc[lambda x: ~np.isin(x['Class'], ['Neuro', 'Per'])]
     .loc[lambda x: ~np.isin(x['Class'], ['Per'])]
     .sort_values(['Class', 'Gene'])
     # .groupby('Class').apply(lambda x: x.sample(frac=.1))
    )

    return cell_genes_merged

def match_cell_genes(cell_genes, pc_version):
    """
    Check which genes are in the PCs
    """
    pc_genes = pc_version.coefs.columns
    cell_types = cell_genes['Class'].unique()
    gene_masks = {}
    for c in cell_types:
        gene_masks[c] = np.isin(pc_genes, cell_genes.query("Class == @c")['Gene'])
    return gene_masks


def shuffle_gene_weights(weights, n=100):
    null_weights = np.repeat(weights.values[:,:,np.newaxis], n, axis=2)
    # null_weights = np.take_along_axis(null_weights, np.random.randn(*null_weights.shape).argsort(axis=0), axis=0)
    for pc in range(3):
        for i in range(n):
            np.random.shuffle(null_weights[:,pc,i])
    
    null_ranks = null_weights.argsort(axis=0).argsort(axis=0)
    nulls = {'weights': null_weights, 'ranks':null_ranks}
    
    return nulls

def compute_cell_scores(weights, nulls, gene_masks, how='mean'):
    """
    Compute cell scores, either mean, or median rank
    """
    true_scores = {}
    null_scores = {}
    
    for name, m in gene_masks.items():
        if how == 'mean':
            true_scores[name] = weights.iloc[m, :].mean(axis=0)
            null_scores[name] = pd.DataFrame(nulls['weights'][m, :, :].mean(axis=0)).T
        elif how == 'median':
            true_ranks = weights.values.argsort(0).argsort(0)
            true_scores[name] = pd.Series(np.median(true_ranks[m, :], axis=0))
            null_scores[name] = pd.DataFrame(np.median(nulls['ranks'][m, :, :], axis=0)).T

    true_scores = pd.concat(true_scores).unstack(1).set_axis(['PC1','PC2','PC3'], axis=1)
    null_scores = pd.concat(null_scores).set_axis(['PC1','PC2','PC3'], axis=1).reset_index(level=0).rename({'level_0':'m'}, axis=1)
    
    return true_scores, null_scores
    
def compute_null_p(true_scores, null_scores, signed=False):
    """
    Compute null p values
    """
    null_p = np.zeros(true_scores.shape)
    for m, name in enumerate(true_scores.index):
        for i in range(3):
            _nulls = null_scores.set_index('m').loc[name].iloc[:,i]
            _score = true_scores.iloc[m, i]
            p = percentileofscore(_nulls, _score)/100
            if signed and p > .5: p = 1-p
            null_p[m, i] = p

    null_p = pd.DataFrame(null_p, 
                          index=true_scores.index, 
                          columns=true_scores.columns)
    
    return null_p