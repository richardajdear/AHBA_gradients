# Functions for enrichment analysis

import numpy as np, pandas as pd
from scipy.stats import percentileofscore
from statsmodels.stats.multitest import multipletests

from processing_helpers import *


### STRING enrichments

def clean_enrichment(file, direction=None, clean_terms=True, filter_001=True):
    """
    Clean STRING enrichment file
    Combine 'Antigen processing' enrichments into one
    """
    enrichment = pd.read_csv(file, delimiter='\t')
    
    if clean_terms:
        # Combine antigen terms
        enrichment.replace({'term description': "Antigen processing.*"}, 'Antigen processing', regex=True, inplace=True)
        enrichment.drop_duplicates('term description', inplace=True)
        # Shorten other names
        replacements = {
            'Negative regulation of cytokine production involved in immune response': 'Regulation of cytokine production',
            'Negative regulation of natural killer cell mediated cytotoxicity': 'Regulation of killer cell mediated cytotoxicity',
            'Protection from natural killer cell mediated cytotoxicity': 'Protection from killer cell mediated cytotoxicity',
            'protein targeting to membrane':'protein targeting',
            'Positive regulation of extrinsic apoptotic signaling pathway': 'Regulation of apoptotic signaling pathway',
            'multicellular organism': '',
            'Establishment of protein localization to': 'Protein localization to',
            'Nuclear-transcribed mrna catabolic process,.*':'mRNA catabolic process'
        }
        enrichment.replace(replacements, inplace=True, regex=True)
        
    
    if direction is not None:
        enrichment = enrichment.loc[lambda x: x['direction'] == direction]
    else:
        enrichment = enrichment.loc[lambda x: x['direction'] != 'both ends']
    
    if filter_001:
        enrichment = enrichment.loc[lambda x: x['false discovery rate'] <= 0.01]
    
    enrichment = (enrichment
                  .assign(neglogFDR = lambda x: -np.log10(x['false discovery rate']))
                  .assign(rank = lambda x: x.groupby('direction')['neglogFDR'].rank(method='first'))
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



def get_cell_genes(which='wen', include=None, subtype=False, combine_layers=False, combine_ex_in=False):
    """
    Read cell genes table
    """

    if which == 'wen':
        path = '../data/wen_cell_genes.csv'
        cell_genes = (
            pd.read_csv(path)
            .set_axis(['Gene','Class'],axis=1)
            .loc[lambda x: np.isin(x['Class'], ['purple','brown','blue'])]
            .replace({'Class':{
                'purple':'Neurons', 'brown':'Synapses', 'blue':'Glia'
            }})
        )
    elif which == 'zeng':
        path = '../data/zeng_layers.csv'
        cell_genes = (
            pd.read_csv(path)
            .set_axis(['Class','Gene'],axis=1)
        )
    elif which == 'jakob':
        path="../data/jakob_cell_genes.csv"
        cell_genes = pd.read_csv(path)

        if include == 'only_lake':
            cell_genes = cell_genes.loc[lambda x: x['Paper'] == 'Lake', :]
        elif include == 'lake_plus_glia':
            cell_genes = cell_genes.loc[lambda x: (x['Paper'] == 'Lake') | (~x['Class'].str.contains('Neuro')), :]
        elif include == 'not_lake':
            cell_genes = cell_genes.loc[lambda x: x['Paper'] != 'Lake', :]

        if subtype:
            cell_genes = (cell_genes
            .assign(Class = lambda x: [t if 'Neuro' in c else c for c,t in zip(x['Class'], x['Type'])])
             )
            
        import re
        if combine_layers:
            cell_genes = (cell_genes
            .assign(Class = lambda x: [re.sub('(a|b|c|d|e)', '', c) if bool(re.search('(Ex|In)', c)) else c for c in x['Class']])
             )
            
        if combine_ex_in:
            cell_genes = (cell_genes
             .replace({'Class':{'Neuro-Ex':'Neuro', 'Neuro-In':'Neuro'}})
            )
        else:
            cell_genes = cell_genes.query("Class != 'Neuro'")
            
        cell_genes = (cell_genes
         # .drop('Class', axis=1).assign(Class = lambda x: x['Type'])
         .melt(id_vars=['Type', 'Paper', 'Cluster', 'Class'], value_name='Gene')
         .loc[lambda x: x['Gene'].notnull(), ['Class', 'Gene']]
         .drop_duplicates()
         # .loc[lambda x: ~np.isin(x['Class'], ['Neuro', 'Per'])]
         .loc[lambda x: ~np.isin(x['Class'], ['Per'])]
         .sort_values(['Class', 'Gene'])
         # .groupby('Class').apply(lambda x: x.sample(frac=.1))
        )

    return cell_genes


def get_cell_genes_weighted():
    """
    x
    """
    
    lake_ex = pd.read_csv("../data/lake_ex.csv")
    lake_in = pd.read_csv("../data/lake_in.csv")
    
    def fisher_norm(y):
        return (np.exp(y) - np.exp(-y))/(np.exp(y) + np.exp(-y))
    
    def clean_lake(df, normalize = True):
        df = (df
         .rename({'cluster':'Class'}, axis=1)
         .melt(id_vars=['Class', 'Gene'], var_name='which', value_name='weight')
         .loc[lambda x: x['Class'] == x['which']].drop('which', axis=1)
         .dropna()
        )
        
        if normalize:
            df = (df
                 .set_index('Gene')
                 # .assign(weight = lambda x: x.groupby('Class').transform(lambda y: (y-np.mean(y))/np.std(y)))
                 .assign(weight = lambda x: x.groupby('Class').transform(lambda y: fisher_norm(y)))
                 .fillna(1)
                 .reset_index()
                 )
        return df
    
    lake_ex = clean_lake(lake_ex, normalize=True)
    lake_in = clean_lake(lake_in, normalize=True)
    
    return pd.concat([lake_ex, lake_in])


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
    
def compute_null_p(true_scores, null_scores, adjust=None, order=None):
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
    
    
    null_p = (null_p
              .reset_index()
              .rename({'level_0':'celltype', 'level_1':'G'}, axis=1)
             )
    
    if order is not None:
        null_p = (null_p
                  .assign(celltype = lambda x: pd.Categorical(x['celltype'],
                                                  ordered=True,
                                                  categories=order))
                  .sort_values('celltype')
             )
    
    return null_p


from sklearn.decomposition import PCA

def make_cell_maps(version, cell_genes, atlas='hcp', normalize = 'std', method='eigengene', order=None, use_weights=True):
    cell_maps = {}
    for celltype in cell_genes.Class.unique():
        cell_genes_ = cell_genes.set_index('Class').loc[celltype]
        if len(cell_genes_)== 1: 
            cell_genes_ = cell_genes_.to_frame().T
        mask = set(cell_genes_.Gene).intersection(version.expression.columns)
        if len(mask) == 0:
            continue
        cell_gene_expression_ = version.expression.loc[:, mask]

        if 'weight' in cell_genes.columns and use_weights:
            cell_gene_weights = cell_genes_.set_index('Gene').loc[mask, 'weight']
            cell_gene_expression_ = cell_gene_expression_ * cell_gene_weights 
        
        # Mean expression
        cell_mean = cell_gene_expression_.mean(axis=1)
        cell_eigengene = PCA(n_components=1).fit_transform(cell_gene_expression_).squeeze()
        cell_eigengene = pd.Series(cell_eigengene, index=cell_gene_expression_.index)
        corr = cell_mean.corr(cell_eigengene)
        # print(f"Class={celltype}: corr(mean expression, eigengene) = {corr}")
        if corr < 0:
            cell_eigengene *= -1
        
        if method=='mean':
            cell_maps[celltype] = cell_mean
        elif method=='eigengene':
            cell_maps[celltype] = cell_eigengene

    cell_maps = pd.concat(cell_maps, axis=1)
    
    if order is not None:
        cell_maps = cell_maps.loc[:,order]
    
    if normalize == 'std':
        cell_maps = cell_maps.apply(lambda x: (x-np.mean(x))/np.std(x))
    elif normalize == 'fisher':
        cell_maps = cell_maps.apply(lambda x: 0.5*(1 + (np.exp(x) - np.exp(-x)) / (np.exp(x) + np.exp(-x)) ))
    
    if atlas == 'hcp':
        cell_maps = cell_maps.join(get_labels_hcp()).set_index('label')
    elif atlas == 'dk':
        cell_maps = cell_maps.join(get_labels_dk()).set_index('label')
        
    return cell_maps


def make_gene_corrs(genes, version, maps):
    genes = set(genes).intersection(version.expression.columns)
    gene_maps = (version.expression.loc[:,genes]
                 .join(get_labels_hcp()).set_index('label'))

    corrs = (pd.concat([gene_maps, maps], axis=1)
             .corr().iloc[:len(genes),len(genes):])
    return corrs
