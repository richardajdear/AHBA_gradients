# Helper functions for analysis

import numpy as np, pandas as pd
from scipy.linalg import orthogonal_procrustes
from itertools import combinations
from processing_helpers import *
from brainsmash.mapgen.base import Base
from statsmodels.formula.api import ols
from sklearn.linear_model import LinearRegression
from scipy import sparse


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


def correlate(a,b):
    """
    Correlate
    """
    n = a.shape[1]
    corr = pd.concat([a,b],axis=1).corr().iloc[:n,n:]
    return corr


def make_corrs_plot(corrs_dict, weights=True):
    """
    Correlate dictionary of version pairs for plotting
    """
    corrs_plot = {}
    corrs_plot['Scores'] = pd.concat({name:pair[1].corr_scores(pair[0]) for name, pair in corrs_dict.items()})
    
    if weights:
        corrs_plot['Weights'] = pd.concat({name:pair[1].corr_weights(pair[0]) for name, pair in corrs_dict.items()})

    corrs_plot = (pd.concat(corrs_plot)
              .stack().to_frame().reset_index()
              .set_axis(['how','version', 'x','y','corr'], axis=1)
             )
    return corrs_plot


def make_var_exp_df(version_dict):
    """
    Get df of variance explained % from dictionary of PCA versions for plotting
    Note: unclear how to calculate this for DME
    """
    dict_var_exp = {name: version.eigenvalues / version.expression.var().sum()
        for name, version in version_dict.items()
    }
    df_var_exp = (pd.DataFrame(dict_var_exp, index=[i+1 for i in range(5)])
        .melt(ignore_index=False).reset_index().set_axis(['PC','version','var'],axis=1)
    )
    return df_var_exp


### Regress out components
def regress_out_components(version, n_components = 1, norm=False):
    data = version.expression
    components = version.scores.iloc[:, 0:n_components]
    lm = LinearRegression().fit(components, data)
    estimated_data = lm.predict(components)
    residuals = data - estimated_data

    if norm:
        residuals = residuals.apply(lambda x: (x-np.mean(x))/np.std(x))
    return residuals


### Get variance explained by regressing out components
def get_var_explained(version):
    original_var = version.expression.var().sum() 
    C1_regressed_var = regress_out_components(version, 1).var().sum()
    C12_regressed_var = regress_out_components(version, 2).var().sum()
    C123_regressed_var = regress_out_components(version, 3).var().sum()

    var_explained_pct = np.array([
        (original_var - C1_regressed_var)/original_var,
        (C1_regressed_var - C12_regressed_var)/original_var,
        (C12_regressed_var - C123_regressed_var)/original_var
    ])
    return var_explained_pct


### Get region-region coexpression matrix
def get_coexp(expression, order='', return_order = False):

    coexp = expression.T.corr()

    if order=='louvain':
        # Sort by louvain modules
        sp = sparse.csr_matrix(coexp+1)
        modules = Louvain(resolution=1.2).fit_transform(sp)
        ord = modules.argsort()
        coexp = coexp.iloc[ord, ord]
    else:
        # Sort by anatomy
        ord = (fetch_hcp()['info'][:180]
              .sort_values(['Lobe', 'cortex'])
              .set_index('id')
              .loc[lambda x: np.isin(x.index, expression.index)]
              .index
            )
        coexp = coexp.loc[ord, ord]
    
    coexp = (coexp
             .set_axis(range(expression.shape[0]))
             .set_axis(range(expression.shape[0]), axis=1)
            )
    
    if return_order:
        return coexp, ord
    else:
        return coexp

### Clustering coefficients
def get_graph_metric(coexp, metric='transitivity', threshold=None, nonnegative=True):
    # sp = sparse.csr_matrix(coexp+1)
    # return Triangles().fit(sp).clustering_coef_
    if nonnegative:
        A = np.where(coexp.values > 0, coexp.values, 0)
    else:
        A = coexp.values

    if threshold is not None:
        # Binarise network
        A = A > np.quantile(A.flatten(), threshold)
    
    G = nx.from_numpy_array(A)

    if metric=='transitivity':
        out = nx.transitivity(G)
    elif metric=='clustering':
        node_coefs = nx.clustering(G, weight='weight')
        out = np.array([*node_coefs.values()]).mean()
    return out


### Smoothness
def get_moran_I(scores,
                pct_include = 25,
                distmat = '../data/parcellations/LeftParcelGeodesicDistmat.txt'):
    # Get distmat matching scores
    inds = [i-1 for i in scores.index]
    distmat=np.loadtxt(distmat)[inds,:][:, inds]

    # # Get indices of region pairs to define distances
    triangle = np.triu_indices(scores.shape[0], k=1)
    distances = distmat[triangle]

    # Filter for distances < percentile
    include = distances < np.percentile(distances, pct_include)
    z_i = scores.values[triangle[0][include]]
    z_j = scores.values[triangle[1][include]]

    # Define weights as 1/distance
    w_ij = 1/distances[include]

    # Compute Moran's I
    I = len(scores) / np.sum(w_ij) * \
        np.sum(w_ij * z_i * z_j) / np.sum(np.square(scores))

    return I

def get_version_variograms(scores, 
                           distmat='../data/LeftParcelGeodesicDistmat.txt',
                           smoothed=False, n_bins=25,
                           return_variograms=True):
    """
    Make variograms for all components in a gradient version
    """
    variograms = {}
    slopes = {}
    for g in ['G1','G2','G3']:
        _scores = scores.loc[:,g]
        variograms[g] = get_variogram(_scores, distmat=distmat, 
                                      smoothed=smoothed, n_bins=n_bins)
        slopes[g] = fit_variogram_slope(variograms[g])
    
    variograms = pd.concat(variograms).reset_index(level=0).rename({'level_0':'G'},axis=1)
    slopes = pd.concat(slopes).unstack()
    if return_variograms:
        return slopes, variograms
    else:
        return slopes

def get_variogram(scores, distmat='../data/LeftParcelGeodesicDistmat.txt', 
                  n_bins=25, pct_include=.25,
                  smoothed=False):
    """
    Make variogram of scores based on distance matrix
    """
    # Filter distmat to match nonmissing scores
    inds = [i-1 for i in scores.index]
    distmat=np.loadtxt(distmat)[inds,:][:, inds]
    
    # Get indices of region pairs
    triangle = np.triu_indices(scores.shape[0], k=1)

    # Select distances from distmat and bin
    x = distmat[triangle]
    bin_counts, bin_edges = np.histogram(x, bins=n_bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
    bin_inds = np.digitize(x, bin_edges[:-1]) # drop rightmost edge to include max value

    # Get squared diff, and bin
    diff_ij = scores.values[triangle[0]] - scores.values[triangle[1]]
    var = np.square(diff_ij)*0.5
    bin_vars = [var[bin_inds==i].mean() for i in np.unique(bin_inds)]

    # Make variogram
    variogram = pd.DataFrame({
        'distance':bin_centers, 
        'variance':bin_vars,
        'counts': bin_counts,
        'include': bin_counts.cumsum()/bin_counts.sum() < pct_include
        })

    # Optionally use smoothing kernel version
    if smoothed:
        mapgen = Base(scores.values, D=distmat, nh=n_bins, pv=100)
        variogram['distance'] = mapgen.h
        variogram['variance'] = mapgen._smvar
    
    return variogram

def fit_variogram_slope(variogram):
    """
    Fit linear model to variogram, ignore last n points
    Return slope and intercept
    """
    variogram_include = variogram.loc[lambda x: x['include']]
    fit = ols(formula='variance ~ distance', data=variogram_include).fit()
    return fit.params


### Enrichment

def sort_genes(df, col):
    out = pd.concat([
        df.sort_values(col, ascending=False).index.to_series().reset_index(drop=True).rename(f'{col}_pos'),
        df.sort_values(col, ascending=True).index.to_series().reset_index(drop=True).rename(f'{col}_neg')
    ], axis=1)
    return out


def sort_genes_PCs(version, i=0, asc=False):
    name = f'PC{i+1}_neg' if asc else f'PC{i+1}_pos'
    return version.coefs.T.iloc[:,:3].sort_values(i, ascending=asc).index.to_series().reset_index(drop=True).rename(name)


def output_PC_gene_ranks(version, name):
    pd.concat([sort_genes_PCs(version, i, asc) for i in range(3) for asc in [False, True]], axis=1).to_csv(f"../outputs/{name}.csv")

def output_PC_region_scores(scores, name, atlas='hcp'):
    if atlas=='hcp':
        labels=get_labels_hcp()[:180]
    elif atlas=='dk':
        labels=get_labels_dk()[:34]
        
    labels.to_frame().join(scores.iloc[:,:3]).set_axis(['label','PC1','PC2','PC3'],axis=1).to_csv(f"../outputs/{name}.csv")

    
    
### Legacy functions from MPhil    
    
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