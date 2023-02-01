# Functions for enrichment analysis

import numpy as np, pandas as pd
from scipy.stats import percentileofscore
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
# import gseapy as gp

from processing_helpers import *


### Functions to compute enrichment from a set of genes, with or without weights


def shuffle_gene_weights(weights, n=100, rank=False):
    """
    Make 'naive' null model by randomizing gene weights / ranks
    """
    null_weights = np.repeat(weights.values[:,:3,np.newaxis], n, axis=2)
    # null_weights = np.take_along_axis(null_weights, np.random.randn(*null_weights.shape).argsort(axis=0), axis=0)
    for g in range(3):
        for i in range(n):
            np.random.shuffle(null_weights[:,g,i])
    
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


def compute_enrichments(weights, null_weights, gene_labels, how='mean', norm=False, posneg=None):
    """
    Compute scores for each gene label, either mean, or median rank
    """
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

    true_enrichments = pd.concat(true_enrichments).unstack(1).set_axis(['G1','G2','G3'], axis=1)
    null_enrichments = pd.concat(null_enrichments).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)

    return true_enrichments, null_enrichments, gene_counts
    
    
def compute_null_p(true_enrichments, null_enrichments, gene_counts=None, adjust='fdr_bh', order=None):
    """
    Compute null p values
    """
    null_pct = np.zeros(true_enrichments.shape)
    for m, label in enumerate(true_enrichments.index):
        for i in range(3):
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
              .rename({'level_0':'label', 'level_1':'G'}, axis=1)
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
        # enrichment.drop_duplicates('term description', inplace=True)
        # Shorten other names
        replacements = {
            'Negative regulation of cytokine production involved in immune response': 'Regulation of cytokine production',
            'Negative regulation of natural killer cell mediated cytotoxicity': 'Regulation of killer cell mediated cytotoxicity',
            'Protection from natural killer cell mediated cytotoxicity': 'Protection from killer cell mediated cytotoxicity',
            'protein targeting to membrane':'protein targeting',
            'Positive regulation of extrinsic apoptotic signaling pathway': 'Regulation of apoptotic signaling pathway',
            'multicellular organism': '',
            'Establishment of protein localization to': 'Protein localization to',
            'Nuclear-transcribed mrna catabolic process,.*':'mRNA catabolic process',
            'Energy derivation by oxidation of organic compounds':'Energy derivation by oxidation',
            'Negative regulation of gene expression, epigenetic':'Regulation of gene expression, epigenetic',
            'Cellular process involved in reproduction in':'Cellular process involved in reproduction',
            'Regulation of dendritic spine development':'Dendritic spine development',
            'Regulation of dendrite development':'Dendrite development',
            'Serotonin receptor signaling pathway':'Serotonin receptor signaling',
            'SRP-dependent cotranslational protein targeting': 'SRP-dependent protein targeting'
        }
        enrichment.replace(replacements, inplace=True, regex=True)
        
        enrichment.drop_duplicates('term description', inplace=True)
    
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


def combine_enrichments(version_, type_, dir_="../outputs/string_data/", include_g1=False, directed=False, top_n=None, **kwargs):
    """
    Combine enrichments from STRING
    """
    directions = {'G1':'top', 'G2':'bottom', 'G3':'bottom'}
    if not directed:
        directions = {k:None for k in directions.keys()}
    
    d = {
        'G2': clean_enrichment(dir_ + version_ + '/g2' + '.enrichment.' + type_ + '.tsv', direction = directions['G2'], **kwargs),
        'G3': clean_enrichment(dir_ + version_ + '/g3' + '.enrichment.' + type_ + '.tsv', direction = directions['G3'], **kwargs)
    }
    
    if include_g1:
        d.update({
            'G1': clean_enrichment(dir_ + version_ + '/g1' + '.enrichment.' + type_ + '.tsv', direction = directions['G1'], **kwargs)
        })
    
    df = (pd.concat(d)
          .reset_index(0).rename({'level_0':'G'}, axis=1)
         )
    
    if directed:
        df['direction'] = 'bottom'
    
    if top_n is not None:
        df = (df
               .assign(rank = lambda x: x.groupby(['G','direction'])['neglogFDR'].rank(method='first',ascending=False))
             .loc[lambda x: x['rank']<=top_n]
             .assign(rank = lambda x: top_n+1-x['rank'])
            )
    
    return df





def get_cell_genes(which='jakob', include=None, subtype=False, combine_layers=False, combine_ex_in=False, add_synapses=False):
    """
    Read cell genes table
    """

    if which == 'wen':
        path = '../data/wen_cell_genes.csv'
        cell_genes = (
            pd.read_csv(path)
            .set_axis(['gene','label'],axis=1)
            .loc[lambda x: np.isin(x['label'], ['purple','brown','blue'])]
            .replace({'label':{
                'purple':'Neurons', 'brown':'Synapses', 'blue':'Glia'
            }})
            .set_index('label').loc[['Neurons','Synapses','Glia'],:].reset_index()
        )
    elif which == 'zeng':
        path = '../data/zeng_layers.csv'
        cell_genes = pd.read_csv(path).set_axis(['label','gene'], axis=1)
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
         .melt(id_vars=['Type', 'Paper', 'Cluster', 'Class'], value_name='gene')
         .loc[lambda x: x['gene'].notnull(), ['Class', 'gene']]
         .rename({'Class':'label'}, axis=1)
         .drop_duplicates()
         .loc[lambda x: ~np.isin(x['label'], ['Per'])]
         .sort_values(['label', 'gene'])
         # .groupby('Class').apply(lambda x: x.sample(frac=.1))
        )

        if add_synapses:
            synapse_genes = get_synapse_genes()
            # cell_genes = cell_genes.loc[lambda x: ~np.isin(x['gene'], synapse_genes['gene'])]
            cell_genes = pd.concat([cell_genes, synapse_genes])
        
    return cell_genes


def get_synapse_genes():
    synapse_genes = pd.read_csv("../data/synaptome_all.csv", usecols=['gene_symbol']).dropna().iloc[:1886]
    synapse_genes = pd.DataFrame({'label':'Synapses', 'gene':synapse_genes['gene_symbol']})
    return synapse_genes

def get_cell_genes_weighted(which=None, normalize=True):
    """
    Read genes table with weights
    """
    
    lake_ex = pd.read_csv("../data/lake_ex.csv")
    lake_in = pd.read_csv("../data/lake_in.csv")
    
    def fisher_norm(y):
        return (np.exp(y) - np.exp(-y))/(np.exp(y) + np.exp(-y))
    
    def clean_lake(df, normalize = True):
        df = (df
         .rename({'cluster':'label', 'Gene':'gene'}, axis=1)
         .melt(id_vars=['label', 'gene'], var_name='which', value_name='weight')
         .loc[lambda x: x['label'] == x['which']].drop('which', axis=1)
         .dropna()
        )
        
        if normalize:
            df = (df
                 .set_index('gene')
                 # .assign(weight = lambda x: x.groupby('label').transform(lambda y: (y-np.mean(y))/np.std(y)))
                 .assign(weight = lambda x: x.groupby('label').transform(lambda y: fisher_norm(y)))
                 .fillna(1)
                 .reset_index()
                 )
        return df
    
    lake_ex = clean_lake(lake_ex, normalize=normalize)
    lake_in = clean_lake(lake_in, normalize=normalize)
    
    if which=='ex':
        return lake_ex
    elif which=='in':
        return lake_in
    else:
        return pd.concat([lake_ex, lake_in])


def get_layer_genes(which='both', add_hse_genes=False):
    he_layers = (pd.read_csv("../data/he_layers.csv")
                .loc[:,['Gene symbol', 'Layer marker in human', 'Log2FC to other layers in human']]
                .set_axis(['gene', 'label', 'log2FC'], axis=1)
                .dropna()
                .loc[lambda x: x['log2FC']>1]
                .sort_values('label')
                .drop(['log2FC'],axis=1)
                )
    
    maynard_data = pd.read_csv("../data/maynard_layers.csv")

    maynard_tstat = (maynard_data
                    .loc[:,['gene', 't_stat_Layer1','t_stat_Layer2','t_stat_Layer3','t_stat_Layer4','t_stat_Layer5','t_stat_Layer6', 't_stat_WM']]
                    .set_index('gene')
                    .set_axis(['L1','L2','L3','L4','L5','L6', 'WM'], axis=1).reset_index()
                    .melt(id_vars='gene', var_name='label', value_name='tstat')
                    .set_index(['gene', 'label'])
                    )

    maynard_layers = (maynard_data
                        .loc[:,['gene', 'fdr_Layer1','fdr_Layer2','fdr_Layer3','fdr_Layer4','fdr_Layer5','fdr_Layer6', 'fdr_WM']]
                        .set_index('gene')
                        .set_axis(['L1','L2','L3','L4','L5','L6', 'WM'], axis=1).reset_index()
                        .melt(id_vars='gene', var_name='label', value_name='fdr')
                        .set_index(['gene', 'label'])
                        .join(maynard_tstat)
                        .loc[lambda x: (x['fdr']<0.05) & (x['tstat']>0)]
                        .reset_index()
                        .sort_values('label')
                        .drop(['tstat', 'fdr'],axis=1)
                        )
    
    he_maynard_layers = pd.concat([
                            he_layers, #.replace({'L6':'L6/WM'}), 
                            maynard_layers.replace({'L6':'L6', 'WM':'L6'})
                        ]).drop_duplicates() 
    

    if which=='he':
        layer_genes = he_layers
    elif which=='maynard':
        layer_genes = maynard_layers
    elif which=='both':
        layer_genes = he_maynard_layers

    hse_genes = ['BEND5',
                'C1QL2',
                'CACNA1E',
                'COL24A1',
                'COL6A1',
                'CRYM',
                'KCNC3',
                'KCNH4',
                'LGALS1',
                'MFGE8',
                'NEFH',
                'PRSS12',
                'SCN3B',
                'SCN4B',
                'SNCG',
                'SV2C',
                'SYT2',
                'TPBG',
                'VAMP1'
                ]
    hse_genes = pd.DataFrame({'gene':hse_genes, 'label':'HSE'})

    if add_hse_genes:
        layer_genes = pd.concat([hse_genes, layer_genes])

    return layer_genes
    # layer_stats = compute_null_p(*compute_enrichments(weights, null_weights, layer_genes))

    # return layer_stats