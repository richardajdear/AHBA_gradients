# Helper functions for analysis

import numpy as np, pandas as pd
# from plots import *
from scipy.linalg import orthogonal_procrustes
from itertools import combinations
from processing_helpers import *

# def sort_genes(df, col):
#     out = pd.concat([
#         df.sort_values(col, ascending=False).index.to_series().reset_index(drop=True).rename(f'{col}_pos'),
#         df.sort_values(col, ascending=True).index.to_series().reset_index(drop=True).rename(f'{col}_neg')
#     ], axis=1)
#     return out


def sort_genes_PCs(version, i=0, asc=False):
    name = f'PC{i+1}_neg' if asc else f'PC{i+1}_pos'
    return version.coefs.T.iloc[:,:3].sort_values(i, ascending=asc).index.to_series().reset_index(drop=True).rename(name)


def output_PC_gene_ranks(version, name):
    pd.concat([sort_genes_PCs(version, i, asc) for i in range(3) for asc in [False, True]], axis=1).to_csv(f"../outputs/{name}.csv")


def test_params(param=None, params_list=None, atlas=None, base=None, **kwargs):
    """
    Compare PCA for different parameter settings
    """
    versions_dict = {}
    
    for p in params_list:
        param_args={param:p}
        expression = get_expression_abagen(atlas=atlas, verbose=0, 
                                           **param_args, **kwargs)
        versions_dict[p] = pcaVersion(expression, message=False)
        print(f'PCA done for param = {p}')

    coef_corrs_dict = {}
    score_corrs_dict = {}
    for p, version in versions_dict.items():
        coef_corrs_dict[p] = base.corr_coefs(version, match=True)[['corr']].T
        score_corrs_dict[p] = base.corr_scores(version, match=True)[['corr']].T

    out = {
        'coef_corrs': pd.concat(coef_corrs_dict).reset_index(level=1, drop=True),
        'score_corrs': pd.concat(score_corrs_dict).reset_index(level=1, drop=True),
        'versions': versions_dict
    }
    return out



def correlate(a,b):
    return pd.concat([a,b],axis=1).corr().iloc[:5,5:]



def get_clusters(expression, coords):
    hc = hier.linkage(expression, 'complete')
    clust_id = hier.cut_tree(hc, n_clusters=2)
    clust_counts = pd.Series(clust_id.squeeze()).value_counts().rename('samples')
    print(clust_counts)
    
    df_cluster_coords = pd.DataFrame(
        data={'cluster':clust_id.squeeze()},
        index=expression.index
    ).join(coords)
    return df_cluster_coords


def compare_matching(version1, version2, annotation, relabel=True):
    df = (
        annotation.set_index('well_id')[[]]
        .join(version1.labels)
        .join(version2.labels, lsuffix='1', rsuffix='2')
        .loc[lambda x: x['label1'].notnull() | x['label2'].notnull()]
        .assign(equal = lambda x: x['label1']==x['label2'],
                match1 = lambda x: x['label1'].notnull(),
                match2 = lambda x: x['label2'].notnull(),
                n=1
               )
        .groupby(['match1', 'match2', 'equal'])
        .agg({'n':'sum'})
        .assign(pct_count = lambda x: x['n']/x['n'].sum())
    )
    if relabel:
        df = (df
         .assign(labels = [
             'Only assigned\nby Abagen',
             'Only assigned\nby Arnatkevic̆iūtė',
             'Assigned to\ndifferent regions',
             'Assigned to\nsame regions'
            ])
        .reset_index(drop=True)
             )
    return df


def disjoint_corrs(triplets_version, how='coefs', match=True, boot=None):
    triplets_list = [list(x) for x in list(combinations(range(6), 3))]
    triplets_names = [''.join(map(str,x)) for x in triplets_list]
    disjoint_triplets = [list(x) for x in combinations(triplets_names,2) if not any(set(list(x[0])).intersection(set(list(x[1]))))]
    
    corrs = {}
    for pair in disjoint_triplets:
#     for pair in all_pairs:
        name = '-'.join(pair)
        pca1 = triplets_version[pair[0]]
        pca2 = triplets_version[pair[1]]
    
        if how=='coefs':
            df_corr = pca1.corr_coefs(pca2, match=match, boot=boot)
        else:
            df_corr = pca1.corr_scores(pca2, match=match, boot=boot)
    
        if match:
            corrs[name] = df_corr['corr'].values
        else:
            corrs[name] = df_corr.pipe(np.diag)
            
    return pd.DataFrame(corrs)



def make_biplot(pca1, pca2, pca1_name='pca1', pca2_name='pca2', title='Rotations Biplot',
                rotate_what='loadings', plot_what='loadings', n=5, dim1=2, dim2=33, 
                atlas='HCP'):
    """
    Make a biplot
    """

    if rotate_what=='loadings':
        A = pca1.loadings
        B = pca2.loadings
    elif rotate_what=='scores':
        A = pca1.U
        B = pca2.U
    
    R1, _ = orthogonal_procrustes(B.iloc[:,:dim1], A.iloc[:,:dim1])
    R2, _ = orthogonal_procrustes(B.iloc[:,:dim2], A.iloc[:,:dim2])
    
    if plot_what=='loadings':
        A = pca1.loadings
        B = pca2.loadings
        text_labels = 'gene_symbol'
    elif plot_what=='scores':
        if atlas=='HCP':
            A = pca1.U.join(get_labels_hcp()).set_index('label').rename_axis('region')
            B = pca2.U.join(get_labels_hcp()).set_index('label').rename_axis('region')
        elif atlas=='DK':
            A = pca1.U.join(get_labels_dk()).set_index('label').rename_axis('region')
            B = pca2.U.join(get_labels_dk()).set_index('label').rename_axis('region')
        text_labels = 'region'
    
    top_A = list(A.reindex(A[1].abs().sort_values(ascending=False).index)[:n].index)
    top_B = list(A.reindex(A[2].abs().sort_values(ascending=False).index)[:n].index)
    top = top_A + top_B
    
    B_rot1 = (B.iloc[:,:dim1] @ R1).set_axis(list(range(1,dim1+1)),axis=1)
    B_rot2 = (B.iloc[:,:dim2] @ R2).set_axis(list(range(1,dim2+1)),axis=1)
    
#     return (R1, R2, top, A, B, B_rot1, B_rot2)
    
    df_biplot = pd.concat({
        f'{pca1_name}': A.loc[lambda x: x.index.isin(top)], 
        f'{pca2_name}': B.loc[lambda x: x.index.isin(top)], 
        f'{pca2_name} rotated in {dim2} PCs': B_rot2.loc[lambda x: x.index.isin(top)],
        f'{pca2_name} rotated in {dim1} PCs': B_rot1.loc[lambda x: x.index.isin(top)],
    })
#     return df_biplot
    
    plot = (df_biplot
        .reset_index()
        .rename(columns={1:'PC1',2:'PC2', 'level_0':'version'})
        .assign(version = lambda x: pd.Categorical(x['version'], categories=x['version'].unique(), ordered=True))
        .pipe(ggplot, aes(x='PC1',y='PC2')) + 
        facet_wrap('~version') +
#           geom_point(size=.2) + 
        geom_text(aes(label=text_labels, color=text_labels),size=8) +
        geom_hline(aes(yintercept=0), size=.2) + geom_vline(aes(xintercept=0), size=.2) +
#         scale_colour_cmap_d() +
        theme_minimal() + coord_fixed() + theme(figure_size=(8,6)) + 
        theme(panel_spacing=.5) + theme(aspect_ratio=1) +
        ggtitle(title)
    )
    return plot
    
        




def corr_rois(version1, version2):
    """
    Get correlations from two versions in same atlas
    """
    scores1 = version1.score_from(version1)
    scores2 = version2.score_from(version2)
    return (
        pd.concat([scores1,scores2],axis=1)
        .corr().iloc[:5,5:]
    )





def compare_matching_by_roi(df):
    aurina_pct = (df.groupby(['label_1']).agg({'equal':['count','sum']}).droplevel(0,axis=1)
                  .assign(
                      n_aur=lambda x: x['count'],
                      diff_aur=lambda x: x['count']-x['sum'],
                      pct_aur=lambda x: 1-x['sum']/x['count'])
                  .drop(['count','sum'],axis=1))
    abagen_pct = (df.groupby(['label_2']).agg({'equal':['count','sum']}).droplevel(0,axis=1)
                  .assign(
                      n_aba=lambda x: x['count'],
                      diff_aba=lambda x: x['count']-x['sum'],
                      pct_aba=lambda x: 1-x['sum']/x['count'])
                  .drop(['count','sum'],axis=1))

    return (pd.DataFrame(index=np.sort(df_dk['label_1'].unique()))
     .join(aurina_pct)
     .join(abagen_pct)
    )




def plot_equivalence(df1, df2, title='Sample matching equivalence, Abagen vs Aurina'):
    df = (
        pd.concat({'DK':df1, 'HCP':df2})
        .reset_index(level=0).rename(columns={'level_0':'atlas'})
    )

    return (ggplot(df) + 
     facet_wrap('~atlas') +
     geom_col(aes(x='reorder(labels, pct_count)', y='pct_count', fill='atlas')) + 
    #  scale_fill_manual([rgb2hex(c) for c in cm.get_cmap('viridis',3).colors][:2][::-1], guide=False) +
     geom_text(aes(x='reorder(labels, pct_count)', y='pct_count', label='n'), ha='right', color='white') + 
     coord_flip() + theme_minimal() + theme(figure_size=(6,2)) +
     xlab('') + ylab('% of samples') +
     ggtitle(title)
    )