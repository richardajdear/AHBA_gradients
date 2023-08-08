# Functions to read and clean data for disorders analysis

import numpy as np, pandas as pd
import mygene
import nibabel as nib
from neuromaps.nulls.spins import parcels_to_vertices, vertices_to_parcels
from processing_helpers import *
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests



def fisher_test(x, target='target', reference='reference', return_table=False):
    """
    Fisher test using dataframe to ensure clean counts
    """
    table = pd.crosstab(x[target], x[reference])
    if return_table:
        print(table)
    odds_ratio, p = fisher_exact(table, 'greater')
    return odds_ratio, p


def test_one_layer(background, disorder_genes, layer_genes,
                   target_layer, reference_disorder, reference_label, 
                   reference_filter = None
                   ):
    """
    Test for enrichment of marker genes of one layer, for a given disorder and label (i.e. GWAS or DEG)
    Optionally add additional filter on the target layer genes (e.g. only C3-positive genes)
    """
    target_layer_genes = layer_genes.loc[lambda x: x['label']==target_layer, 'gene']
    reference_disorder_genes = (disorder_genes
                                .loc[lambda x: (x['disorder']==reference_disorder) & 
                                               (x['label']==reference_label), 'gene'])
    
    df = (pd.DataFrame(index=background)
        .assign(is_reference = lambda x: np.isin(x.index, reference_disorder_genes))
        .assign(is_target = lambda x: np.isin(x.index, target_layer_genes))
    )

    if reference_filter is not None:
        df = df.assign(is_reference = lambda x: x['is_reference'] & np.isin(x.index, reference_filter))

    n_in_reference = sum(df['is_reference'])
    n_target_in_reference = sum(df['is_target'] & df['is_reference'])
    odds_ratio, p = df.pipe(fisher_test, target='is_target', reference='is_reference', return_table=False)
    return n_in_reference, n_target_in_reference, odds_ratio, p


def test_layers_all_combinations(background, disorder_genes, layer_genes, reference_filter):
    """
    Test layer enrichments of all combinations of GWAS and DEG genes for multiple disorders
    Optionally add additional filter on the target layer genes (e.g. only C3-positive genes)
    """
    enrichments_list = []
    for disorder in disorder_genes['disorder'].unique():
        unique_labels = disorder_genes.loc[lambda x: x['disorder']==disorder, 'label'].unique()
        for label in unique_labels:
            for layer in layer_genes['label'].unique():
                # print(f"Testing {label}-{disorder}-{layer}...")
                n_in_disorder, n_layer_in_disorder, odds_ratio, p = test_one_layer(
                                                                background = background, 
                                                                disorder_genes = disorder_genes,
                                                                layer_genes = layer_genes,
                                                                target_layer = layer, 
                                                                reference_disorder = disorder, 
                                                                reference_label = label, 
                                                                reference_filter = reference_filter)
                out = pd.Series({
                        'label': label,
                        'disorder': disorder,
                        'layer': layer,
                        'n': n_layer_in_disorder,
                        'N': n_in_disorder,
                        'pct': n_layer_in_disorder / n_in_disorder,
                        'odds_ratio': odds_ratio,
                        'p': p
                        })
                enrichments_list.append(out)
            print(f"Done {label}-{disorder}.")

    stats = (pd.DataFrame(enrichments_list)
          .assign(q = lambda x: multipletests(x['p'], method='fdr_bh')[1])
          .assign(sig = lambda x: x['q']<0.05)
    )
    return stats



def make_component_quantiles(weights, q=10, labels=None, levels=None, na_value=np.NaN):
    """
    Split weights by component quantiles, then match to labels
    """
    weights_with_quantiles = (weights
            .melt(ignore_index=False, var_name='C', value_name='C_score')
            .reset_index().rename({'index':'gene'}, axis=1)
            .assign(C_quantile = lambda x: x.groupby('C').apply(lambda y: pd.qcut(y['C_score'], q=q, labels=range(1,q+1))).reset_index(0, drop=True))
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
                # .assign(rank_in_quantile = lambda x: x.groupby(['C','C_quantile']).cumcount()+1)
                # .sort_values(['C_quantile','rank_in_quantile'])
            )


def make_label_quantiles_by_component(weights, labels, q=10):
    """
    Split labels into quantiles by component weights
    """
    label_percentiles = (labels.loc[:, ['label','gene']]
        .join(weights.melt(ignore_index=False, var_name='C',value_name='C_score'), on='gene')
        .dropna()
        .reset_index(drop=True)
        .assign(C_quantile = lambda x: x.groupby(['C','label'])['C_score'].apply(lambda y: pd.qcut(y, q=q, labels=range(1,q+1
        ))))
    )
    return label_percentiles



# def aggregate_disorder_gene_quantiles_over_layers(quantiles, layers, disorder='SCZ'):
#     quantile_layers = (quantiles
#         .join(layers.set_index('gene')['label'].rename('layer'), on='gene')
#         .fillna({'layer':'no_layer'})
#         .groupby(['label','C','C_quantile','layer']).count()
#         .assign(pct = lambda x: x['gene']/x.groupby(['label','C','C_quantile'])['gene'].sum())
#         .reset_index()
#         .loc[lambda x: x['C']=='C3']
#         .loc[lambda x: x['layer']!='no_layer']
#         .replace({'label':{'DEG':f'{disorder} RNAseq', 'GWAS':f'{disorder} GWAS'}})
#         )
#     return quantile_layers


# def aggregate_disorder_gene_quantiles_over_development(quantiles, curves, disorder='SCZ'):
#     quantile_curves_disorders = (quantiles
#      .join(curves.set_index('gene'), on='gene').dropna()
#      .groupby(['label','C','C_quantile','age_log10'])
#      .agg({'pred':'mean'})
#      .reset_index()
#      .loc[lambda x: x['C']=='C3']
#      .replace({'label':{'DEG':f'{disorder} RNAseq', 'GWAS':f'{disorder} GWAS'}})
#     )
#     return quantile_curves_disorders

