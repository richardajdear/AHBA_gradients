# Functions for enrichment analysis

import numpy as np, pandas as pd
from scipy.stats import percentileofscore
from statsmodels.stats.multitest import multipletests


### STRING enrichments

def clean_enrichment(file, direction=None, combine_antigen=True):
    """
    Clean STRING enrichment file
    Combine 'Antigen processing' enrichments into one
    """
    enrichment = pd.read_csv(file, delimiter='\t')
    
    if combine_antigen:
        enrichment.replace({'term description': "Antigen processing.*"}, 'Antigen processing', regex=True, inplace=True)
        enrichment.drop_duplicates('term description', inplace=True)
    
    if direction is not None:
        enrichment = enrichment.loc[lambda x: x['direction'] == direction]
    else:
        enrichment = enrichment.loc[lambda x: x['direction'] != 'both ends']
    
    enrichment = (enrichment
                  .assign(neglogFDR = lambda x: -np.log10(x['false discovery rate']))
                  .assign(rank = lambda x: len(x['neglogFDR']) - x['neglogFDR'].rank())
                  .assign(FDR = lambda x: x['false discovery rate'],
                          enrichment = lambda x: x['enrichment score'],
                          n_genes = lambda x: x['genes mapped'],
                          description = lambda x: x['term description']
                         )
                  .loc[:, ['description', 'n_genes', 'direction', 'enrichment', 'FDR', 'neglogFDR', 'rank']]
                 )
    return enrichment


def combine_enrichments(version_, type_, dir_="../outputs/string_data/", include_g1=False, directed=False):
    """
    Combine enrichments
    """
    directions = {'G1':'top', 'G2':'bottom', 'G3':'bottom'}
    if not directed:
        directions = {k:None for k in directions.keys()}
    
    d = {
        'G2': clean_enrichment(dir_ + version_ + '/g2' + '.enrichment.' + type_ + '.tsv', direction = directions['G2']),
        'G3': clean_enrichment(dir_ + version_ + '/g3' + '.enrichment.' + type_ + '.tsv', direction = directions['G3'])
    }
    
    if include_g1:
        d.update({
            'G1': clean_enrichment(dir_ + version_ + '/g1' + '.enrichment.' + type_ + '.tsv', direction = directions['G1'])
        })
    
    df = (pd.concat(d)
          .reset_index(0).rename({'level_0':'G'}, axis=1)
         )
    
    if directed:
        df['direction'] = 'bottom'
    
    return df



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

def match_cell_genes(cell_genes, weights):
    """
    Check which genes are in the PCs
    """
    genes = weights.index
    cell_types = cell_genes['Class'].unique()
    gene_masks = {}
    for c in cell_types:
        gene_masks[c] = np.isin(genes, cell_genes.query("Class == @c")['Gene'])
    return gene_masks


def shuffle_gene_weights(weights, n=100, rank=False):
    null_weights = np.repeat(weights.values[:,:,np.newaxis], n, axis=2)
    # null_weights = np.take_along_axis(null_weights, np.random.randn(*null_weights.shape).argsort(axis=0), axis=0)
    for g in range(3):
        for i in range(n):
            np.random.shuffle(null_weights[:,g,i])
    
    nulls = {'weights': null_weights}
    
    if rank:
        null_ranks = null_weights.argsort(axis=0).argsort(axis=0)
        nulls.update({'ranks':null_ranks})
    
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

    true_scores = pd.concat(true_scores).unstack(1).set_axis(['G1','G2','G3'], axis=1)
    null_scores = pd.concat(null_scores).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'m'}, axis=1)
    
    return true_scores, null_scores
    
def compute_null_p(true_scores, null_scores, adjust=None):
    """
    Compute null p values
    """
    null_pct = np.zeros(true_scores.shape)
    for m, name in enumerate(true_scores.index):
        for i in range(3):
            _nulls = null_scores.set_index('m').loc[name].iloc[:,i]
            _score = true_scores.iloc[m, i]
            pct = percentileofscore(_nulls, _score)/100
            null_pct[m, i] = pct
            
    true_mean = true_scores.stack().rename('true_mean')

    null_mean = (null_scores
                 .groupby('m').agg(['mean','std'])
                 .stack(0)
                 .rename_axis([None,None])
                 .set_axis(['null_mean', 'null_std'], axis=1)
                )
            
    null_p = (pd.DataFrame(null_pct, 
                           index=true_scores.index, 
                           columns=true_scores.columns)
              .stack().rename('pct').to_frame()
              .join(true_mean)
              .join(null_mean)
              .assign(z = lambda x: (x['true_mean'] - x['null_mean'])/x['null_std'])
              .assign(pos = lambda x: [pct > 0.5 for pct in x['pct']])
              .assign(p = lambda x: [(1-pct)*2 if pct>0.5 else pct*2 for pct in x['pct']]) # x2 for bidirectional
              
             )
    
    if adjust is not None:
        null_p = (null_p
             .assign(q = lambda x: multipletests(x['p'], method=adjust)[1])
             .assign(sig = lambda x: x['q'] < .05)
             # .assign(q_abs = lambda x: [1-q if pos else q for pos, q in zip(x['pos'], x['q'])])
            )
    
    order = ['Neuro', 'Oligo', 'Astro', 'Micro', 'Endo', 'OPC']
    null_p = (null_p
              .loc[order]
              .reset_index()
              .rename({'level_0':'celltype', 'level_1':'G'}, axis=1)
              .assign(celltype = lambda x: pd.Categorical(x['celltype'],
                                                  ordered=True,
                                                  categories=order))
             )
    
    return null_p
