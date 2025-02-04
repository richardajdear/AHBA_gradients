# Code for comparing PC maps with MRI maps
import numpy as np, pandas as pd
from processing import *
from analysis_helpers import *

# from netneurotools import freesurfer as nnsurf
from neuromaps.images import annot_to_gifti
from neuromaps.nulls.spins import parcels_to_vertices, vertices_to_parcels
from neuromaps.stats import compare_images
# from scipy.stats import percentileofscore
from statsmodels.stats.multitest import multipletests
from nibabel.freesurfer.io import read_annot
# from brainsmash.mapgen.base import Base
from sklearn.decomposition import PCA
# from nilearn.input_data import NiftiLabelsMasker
import bct


def make_correlation_matrix(maps, gamma=1, seed=1, method='pearson', cluster_reorder_dict={}):
    rename_dict = {
        'T1T2': 'T1/T2',
        'paquola_FC': 'FC',
        'glasser_GI': 'AG',
        'MEG_theta': 'MEGθ',
        'dMT': 'ΔMT',
        'd_intramod60': 'ΔIntraFC',
        'MIND_CT': 'MIND',
        'hill_evo': 'Evo Exp',
        'hill_dev': 'Dev Exp',
        'externopyramidisation': 'Ext',
        'thickness': 'CT',
        'dCT': 'ΔCT',
        'PC1_neurosynth': 'NS PC1',
        'G1_fMRI': 'fMRI G1',
        'allom': 'Allom',
        'CBF': 'CBF',
    }

    corrmat = maps.rename(rename_dict, axis=1).corr(method=method)
    C1_ranks = corrmat.abs()['C1'].rank(ascending=False)
    C2_ranks = corrmat.abs()['C2'].rank(ascending=True)
    C3_ranks = corrmat.abs()['C3'].rank(ascending=True)
    combo_ranks = C1_ranks + C3_ranks

    corrmat = corrmat.iloc[combo_ranks.argsort(),:].iloc[:,combo_ranks.argsort()]

    clusters = bct.community_louvain(corrmat.abs().values, gamma=gamma, seed=seed)[0]
    print('n clusters', len(np.unique(clusters)))

    clusters = pd.Series(clusters, index=corrmat.index).replace(cluster_reorder_dict)

    corrmat = corrmat.iloc[np.argsort(clusters), :].iloc[:, np.argsort(clusters)]

    corrmat_for_plot = (corrmat
            .rename_axis('x')
            .melt(ignore_index=False, var_name='y', value_name='r')
            .reset_index()
            .assign(cluster = lambda x: x['x'].map(clusters))
    )
    return corrmat_for_plot




def hcp_to_dk(scores, hcp_img = None, dk_img = None, as_df=True):
    """
    Project HCP scores to DK
    """
    if hcp_img is None:
        hcp_img_path = "../data/parcellations/lh.HCPMMP1.annot"
        hcp_img = annot_to_gifti(hcp_img_path)
    if dk_img is None:
        dk_img_path = "../data/parcellations/lh.aparc.annot"
        dk_img = annot_to_gifti(dk_img_path)

    scores_dk = np.zeros((34,3))
    for i in range(3):
        # Re-index gradient null values with NA
        if as_df:
            c_hcp = scores[i].reindex(range(1,181)).values
        else:
            c_hcp = scores[:,i]
        # Use HCP parcellation image to project HCP data to fsaverage
        c_fsaverage = parcels_to_vertices(c_hcp, hcp_img)
        # Use DK parcellation image to project fsaverage data into DK
        c_dk = vertices_to_parcels(c_fsaverage, dk_img)
        # Add to outputs
        scores_dk[:,i] = c_dk

    # Convert to dataframe
    if as_df:
        scores_dk = pd.DataFrame.from_records(scores_dk, index=list(range(1,35)))

    return scores_dk



def dk_to_hcp(maps,
                hcp_img_path = "../data/parcellations/lh.HCPMMP1.annot",
                dk_img_path = "../data/parcellations/lh.aparc.annot"
                ):
    """
    Project DK scores to HCP
    """
    hcp_img = annot_to_gifti(hcp_img_path)
    dk_img = annot_to_gifti(dk_img_path)
    n_maps = maps.shape[1]
    
    maps_hcp = np.zeros((180,n_maps))
    for i in range(n_maps):
        g_dk = maps.iloc[:,i].values
        # Use DK parcellation image to project DK data to fsaverage
        g_fsaverage = parcels_to_vertices(g_dk, dk_img)
        # Use HCP parcellation image to project fsaverage data into HCP
        g_hcp = vertices_to_parcels(g_fsaverage, hcp_img)
        # Add to outputs
        maps_hcp[:,i] = g_hcp
    
    # Convert to dataframe
    maps_hcp = pd.DataFrame.from_records(maps_hcp, 
                                         index=list(range(1,181)),
                                         columns=maps.columns
                                         )
    
    return maps_hcp



def maps_pca(maps, short_names=None, n_components=3):
    """
    Run PCA on maps
    """
    cols = [f'MRI PC{i+1}' for i in range(n_components)]
    pca = PCA(n_components=n_components).fit(maps)
    pca_scores = (pd.DataFrame(pca.transform(maps), 
                               index=maps.index, 
                               columns=cols
                              )
                  .apply(lambda x: (x-np.mean(x))/np.std(x))
                 )
    pca_weights = pd.DataFrame(pca.components_.T, 
                               columns=cols, 
                               index=maps.columns)
    
    pca_var = pca.explained_variance_ratio_
    
    if short_names is not None:
        pca_weights = pca_weights.assign(map_name = lambda x: x.index.map(short_names))
    
    return pca_scores, pca_var, pca_weights


def maps_pca_boot(version, maps, n_boot=10, n_sample=5, n_components=5, replace=False):
    """
    Test maps PCA over random permutations
    """
    pca_corrs_boot = {}
    for i in range(n_boot):
        maps_boot = maps.sample(n=n_sample, replace=replace, axis=1)
        pca_scores, pca_var, pca_weights = maps_pca(maps_boot, n_components=n_components)
        version_scores = version.clean_scores(n_components=n_components).set_index('label')
        corrs = pca_scores.join(scores).corr().iloc[:n_components,n_components:].abs()
        matches = version.match_components(corrs, n_components=n_components)
        index_ = [f'MRI PC{i+1}' for i in range(n_components)]
        var = pd.Series(pca_var, index=index_, name='var')

        MRI_PC = [ix[0] for ix in matches[0]]
        C = [ix[1] for ix in matches[0]]
        pca_corrs_boot[i] = (pd.DataFrame({
                                'MRI_PC':MRI_PC, 
                                'C':C, 
                                'boot':i,
                                'r':matches[1]
                                })
                        .set_index('MRI_PC').join(var)
                        .reset_index()
                        )

    pca_corrs_boot = pd.concat(pca_corrs_boot).reset_index(drop=True)
    return pca_corrs_boot


# Atlas regions enrichment
def compute_region_means(scores, null_scores, regions_series, method='median', reindex=True):
    if reindex:
        scores = scores.reindex(range(1,181))

    masks = {}
    # regions_series = regions.drop('label',axis=1).squeeze().dropna()
    for region in regions_series.dropna().unique():
        region_ids = regions_series.loc[lambda x: x==region].index
        masks[region] = np.isin(scores.index, region_ids)

    scores_array = scores.iloc[:,:3].values
    
    true = {}
    null = {}
    
    for region, mask in masks.items():
        if method=='mean':
            true[region] = pd.Series(np.nanmean(scores_array[mask, :], axis=0))
            null[region] = pd.DataFrame(np.nanmean(null_scores[mask, :, :], axis=0)).T
        elif method=='median':
            true[region] = pd.Series(np.nanmedian(scores_array[mask, :], axis=0))
            null[region] = pd.DataFrame(np.nanmedian(null_scores[mask, :, :], axis=0)).T                
    
    true = pd.concat(true).unstack(1).set_axis(['C1','C2','C3'], axis=1)
    null = pd.concat(null).set_axis(['C1','C2','C3'], axis=1) \
                 .reset_index(level=0).rename({'level_0':'label'}, axis=1)

    return true, null






# LEGACY
# def generate_spins_from_gradients(scores, n=10,
#                            outfile='../outputs/permutations/spin_gradients_1000.npy',
#                            atlas='hcp'):
    
#     if atlas == 'dk':
#         atlas = 'aparc'
#         _,_,rh_names = read_annot(f"../data/parcellations/rh.{atlas}.annot")
#         drop = rh_names + ['lh_corpuscallosum', 'lh_unknown']
#     else: 
#         atlas = 'HCPMMP1'
#         _,_,rh_names = read_annot("../data/parcellations/rh.HCPMMP1.annot")
#         # Find regions missing in scores
#         missing_rois = list(set(get_labels_hcp()[:180]).difference(scores['label']))
#         # Fix 7Pl to match HCP codes
#         missing_rois = ['7PL' if x=='7Pl' else x for x in missing_rois]
#         # Drop missing regions
#         rh_names = rh_names + [np.bytes_('L_' + roi + '_ROI') for roi in missing_rois]
#         drop = rh_names
    
#     spin_pcs = nnsurf.spin_data(
#         data = np.array(scores.set_index('label')),
#         drop = drop,
#         version = "fsaverage",
#         lhannot = f"../data/parcellations/lh.{atlas}.annot",
#         rhannot = f"../data/parcellations/rh.{atlas}.annot",
#         n_rotate = n,
#     )
#     np.save(outfile, spin_pcs)

###
# def get_corrs(scores, maps, method='pearson', atlas='hcp', n_components=3):
    
#     # labels = get_labels_dk() if atlas=='dk' else get_labels_hcp()[:180]
        
#     corrs = (
#         scores
#         .set_index('label')
#         .join(maps)
#         .corr(method=method).iloc[n_components:,:n_components]
#         .set_axis([f'C{i+1}' for i in range(n_components)], axis=1)
#     )
#     return corrs


# Legacy method

# def corr_nulls_from_grads(null_grads, scores, maps, method='pearson', pool=False, pool_frac=.1):
#     """
#     Get correlations with maps from gradient score nulls
#     """            
#     # Filter maps for regions in gradients
#     maps_filter = maps.set_axis(range(1, maps.shape[0]+1)).loc[scores.index, :]

#     null_corrs = {}
#     for m in range(null_grads.shape[1]):
#         # Optionally pool maps together
#         if pool:
#             # Set dimensions for reshaping
#             # Downsample nulls because of pooling
#             m = null_grads.shape[1]
#             n = int(null_grads.shape[2] * pool_frac)
#             null_grads_downsample = null_grads[:,:,:n]
#             nulls = null_grads_downsample.reshape(-1, m*n)
#         else:
#             nulls = null_grads[:,m,:]
        
#         # Concat scores and nulls
#         # df_concat = pd.concat([maps.set_axis(maps_index), nulls], axis=1)
#         # Optionally correlate with spearman, otherwise pearson
#         if method == 'spearman':
#             df_corr = np_pearson_corr(maps_filter.rank(), nulls.rank())
#         else:
#             df_corr = np_pearson_corr(maps_filter, nulls)
#         # Cleanup and stack maps into vector
#         null_corrs[m] = df_corr.stack().droplevel(1)
        
#     # Concat columns (each gradient)
#     null_corrs = (
#         pd.concat(null_corrs, axis=1).reset_index(level=0)
#         .set_axis(['map','G1','G2','G3'], axis=1)
#     )
#     return null_corrs
    

# def corr_nulls_from_maps(null_maps, scores, maps, method='pearson', pool=False, pool_frac=.1):
#     """
#     Get correlations with scores from null maps
#     """
    
#     # Drpp 'label' field if it exists
#     if 'label' in scores.columns:
#         scores = scores.drop('label', axis=1)
        
#     if pool:
#         # Set dimensions for reshaping
#         # Downsample nulls because of pooling
#         m = null_maps.shape[1]
#         n = int(null_maps.shape[2] * pool_frac)
#         null_maps_downsample = null_maps[:,:,:n]
#         null_maps_pool = null_maps_downsample.reshape(-1, m*n)
    
#     null_corrs = {}
#     # Define index as 1 to n
#     index = [i+1 for i in range(null_maps.shape[0])]
#     for m, mapname in enumerate(maps.columns):
#         # Optionally pool maps together
#         if pool:
#             nulls = pd.DataFrame(null_maps_pool, index=index)
#         else:
#             nulls = pd.DataFrame(null_maps[:,m,:], index=index)
        
#         # Concat scores and nulls, matching on index directly
#         df_concat = pd.concat([scores, nulls], axis=1)
#         # Optionally correlate with spearman, otherwise pearson
#         if method == 'spearman':
#             # .rank().corr() is faster than .corr(method='spearman') because of null checks
#             df_corr = df_concat.rank().corr()
#         else:
#             df_corr = df_concat.corr()
#         # Cleanup
#         null_corrs[mapname] = df_corr.iloc[3:,:3]
        
#     # Concat rows (each map)
#     null_corrs = (
#         pd.concat(null_corrs).reset_index(level=0)
#         .set_axis(['map','G1','G2','G3'],axis=1)
#     )
#     return null_corrs


# def get_null_p(true_corrs, null_corrs, adjust='fdr_bh', order=None):
#     """
#     Get p values
#     """
#     # For each map, for each set of nulls, compute the percentile of the true correlations
#     null_pct = np.zeros(true_corrs.shape)
#     for m, _map in enumerate(true_corrs.index):
#         for i in range(3):
#             _null_corrs = null_corrs.set_index('map').loc[_map].iloc[:,i]
#             _corr = true_corrs.iloc[m,i]
#             pct = percentileofscore(_null_corrs, _corr)/100
#             null_pct[m,i] = pct

#     true_mean = true_corrs.stack().rename('true_mean')

#     null_mean = (null_corrs
#                  .groupby('map').agg(['mean','std'])
#                  .stack(0)
#                  .rename_axis([None,None])
#                  .set_axis(['null_mean', 'null_std'], axis=1)
#                 )
            
#     null_p = (pd.DataFrame(null_pct, index=true_corrs.index, columns=true_corrs.columns)
#               .stack().rename('pct').to_frame()
#               .join(true_mean)
#               .join(null_mean)
#               .assign(z = lambda x: (x['true_mean'] - x['null_mean'])/x['null_std'])
#               .assign(pos = lambda x: [pct > 0.5 for pct in x['pct']])
#               .assign(p = lambda x: [(1-pct)*2 if pct>0.5 else pct*2 for pct in x['pct']]) # x2 for bidirectional
#               .assign(sig = lambda x: x['p'] < .05)
#              )
    
#     if adjust is not None:
#         null_p = (null_p
#              .assign(q = lambda x: multipletests(x['p'], method=adjust)[1])
#              .assign(sig = lambda x: x['q'] < .05)
#              # .assign(q_abs = lambda x: [1-q if pos else q for pos, q in zip(x['pos'], x['q'])])
#             )
#     else:
#         null_p = (null_p
#              .assign(q = lambda x: x['p'])
#              .assign(sig = lambda x: x['p'] < .05)
#                  )
    
#     null_p = (null_p
#               .reset_index()
#               .rename({'level_0':'map', 'level_1':'G'}, axis=1)
#               # .assign(map = lambda x: pd.Categorical(x['map'], ordered=True, categories=order))
#              )
    
#     # Fix order of maps
#     if order is None:
#         order = true_corrs.index
#     null_p = (null_p
#               .assign(map = lambda x: pd.Categorical(x['map'], ordered=True, categories=order))
#               .sort_values('map')
#          )

#     return null_p


    

# def print_corrs_sig(corrs, null_p):
#     map_corrs_sig = (corrs
#      .loc[null_p.index]
#      .round(2).astype('string')
#      .where(null_p > .05, other = lambda x: x+' *')
#      .where(null_p > .01, other = lambda x: x+'*')
#      .where(null_p > .001, other = lambda x: x+'*')
#     )
#     # map_corrs_sig.to_csv("../outputs/map_corrs_sig.csv")
#     return map_corrs_sig


# def get_inflammation_data(img="../data/inflammation/Activation_proportion.nii.gz"):
#     hcp_mask = "../data/parcellations/HCP-MMP_1mm.nii.gz"
#     mask = NiftiLabelsMasker(hcp_mask, resampling_target='data')
#     inflammation_img = img
#     inflammation_hcp = mask.fit_transform(inflammation_img).squeeze()

#     inflammation_hcp = (
#         pd.Series(inflammation_hcp)
#         .rename('activation')
#         .to_frame()
#         .set_axis([i+1 for i in range(360)])
#         .join(get_labels_hcp())
#         .set_index('label')
#         # .rename({'label':'region'}, axis=1)
#         # set_index('region')
#         # .assign(hemi=lambda x:['left' if id<=180 else 'right' for id in x.index])
#                )
#     return inflammation_hcp