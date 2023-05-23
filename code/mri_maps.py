# Code for comparing PC maps with MRI maps
import numpy as np, pandas as pd
from processing_helpers import *
from analysis_helpers import *

from netneurotools import freesurfer as nnsurf
from neuromaps.stats import compare_images
# from scipy.stats import percentileofscore
from statsmodels.stats.multitest import multipletests
from nibabel.freesurfer.io import read_annot
# from brainsmash.mapgen.base import Base
from sklearn.decomposition import PCA
# from nilearn.input_data import NiftiLabelsMasker
import bct



def get_maps(data_dir="../data/cortical_maps.csv"):
    names = {
        'paquola_FC': 'Functional Hubs (FC)',
        'T1T2': 'Myelination (T1w/T2w)',
        'MEG_theta': 'MEG Theta (MEGθ)',
        'glasser_GI': 'Aerobic Glycolysis (AG)',
        'd_intramod60': 'Age 14-26 change in Short-Range\nIntramodular Connectivity (ΔShIm)',
        'dMT': 'Age 14-26 change in \n Myelination (ΔMT)',
    }

    maps = (pd.read_csv(data_dir)
            .apply(lambda x: (x-np.mean(x))/np.std(x))
            .set_index(get_labels_hcp()[:180])
    )

    return maps

def get_mesulam_ve_yeo(data_dir="../data/mesulam_ve_yeo.csv"):
    short_names_dict = {
        'Mesulam_names': {
            'Heteromodal':'Het.',
            'Unimodal':'Uni.',
            'Paralimbic':'Par.',
            'Idiotypic':'Idi.'
        },
        'Yeo_names': {
            'Visual':'VIS',
            'Somatomotor':'SMN',
            'Dorsal attention':'DAN',
            'Central attention':'CAN',
            'Default mode':'DMN',
            'Frontoparietal':'FPN',
            'Limbic':'LIM'
        },
    }

    mesulam_ve_yeo = (pd.read_csv(data_dir, index_col=0)
                      .drop(['x','y','z'], axis=1)
                      .replace({'label':{'7PL':'7Pl'}})
                      .replace({'Mesulam':{4:1,2:4,3:2}})
                    #   .replace({'Yeo_names':{'Central attention':'Central atten.','Dorsal attention':'Dorsal atten.'}})
                      .replace(short_names_dict)
                      .set_index('label'))
    return mesulam_ve_yeo


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
    G1_ranks = corrmat.abs()['G1'].rank(ascending=False)
    G2_ranks = corrmat.abs()['G2'].rank(ascending=True)
    G3_ranks = corrmat.abs()['G3'].rank(ascending=True)
    combo_ranks = G1_ranks + G3_ranks

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



# def get_mri_maps_OLD(data_dir="../data/stat_maps_HCP_forRichard.csv", filter=True):
#     maps = (
#         pd.read_csv(data_dir)
#         .apply(lambda x: (x-np.mean(x))/np.std(x))
#         .sort_index(axis=1)
#         .set_index(get_labels_hcp()[:180])
#         # .rename_axis('region')#.reset_index()
#     )
    
#     selected_maps = [
#         'T1T2',
#         'thickness',
#         'glasser_CMRO2',
#         'glasser_CMRGlu',
#         'G1_fMRI',
#         'PC1_neurosynth',
#         'externopyramidisation',
#         'glasser_GI',
#         'hill.evo_remapped',    
#         'hill.dev_remapped',
#         'glasser_CBF',
#         'allom',
        
#         # 'glasser_CBV',
#         # 'asl',
        
#         # 'x':'X axis',
#         # 'y':'Y axis',
#         # 'z':'Z axis'
#     ]
    
#     if filter:
#         maps = maps.loc[:, selected_maps]
        
#     return maps


def get_ctmt_from_dk(project_to_hcp=True):
    ctmt = (
        pd.read_csv("../data/whitakervertes2016_complete.csv", index_col=0)
        .query('hemi=="l"')
        .set_index('label')
        .loc[:, ['CT','CT_delta','MT','MT_delta']]#,'PLS2']]
    )
    
    if project_to_hcp:
        ctmt = ctmt.pipe(dk_to_hcp).set_axis(get_labels_hcp()[:180])
            
    return ctmt


def get_disorder_maps(data_dir="../data/lifespan_dx_DKatlas.csv"):
    maps = (
        pd.read_csv(data_dir, index_col=0)
        # .apply(lambda x: (x-np.mean(x))/np.std(x))
        .sort_index(axis=1)
        .rename_axis('label')#.reset_index()
    )
    maps = maps.set_index(maps.index.str.replace(" ", ""))
    maps = maps.set_index('lh_' + maps.index)
 
    return maps

def get_brainchart_maps(data="../data/Peaks_Table_2_2.csv"):
    maps = (
        pd.read_csv(data, index_col=0)
        .drop('Nboot',axis=1)
        .set_index('feat')
        .rename_axis('label')
    )
    maps = maps.set_index('lh_' + maps.index)
 
    return maps    


def get_corrs(scores, maps, method='pearson', atlas='hcp', n_components=3):
    
    # labels = get_labels_dk() if atlas=='dk' else get_labels_hcp()[:180]
        
    corrs = (
        scores
        .set_index('label')
        .join(maps)
        .corr(method=method).iloc[n_components:,:n_components]
        .set_axis([f'G{i+1}' for i in range(n_components)], axis=1)
    )
    return corrs



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



def generate_shuffles(maps, n=1000,
                      outfile='../outputs/shuffle_maps_1000.npy'):
    
    if 'label' in maps.columns:
        maps = maps.drop('label', axis=1)
    
    shuffle_maps = np.repeat(maps.values[:,:,np.newaxis], 1000, axis=2)
    for i in range(shuffle_maps.shape[2]):
        np.random.shuffle(shuffle_maps[:,:,i])

    np.save(outfile, shuffle_maps)


def generate_spins(maps, n=10, 
                   outfile='../outputs/spin_maps_1000.npy',
                   atlas='hcp'):
    """
    Generate spins from a set of maps and save to file
    """
    
    if atlas == 'dk':
        atlas = 'aparc'
        _,_,rh_names = read_annot(f"../data/parcellations/rh.{atlas}.annot")
        drop = rh_names + ['lh_corpuscallosum', 'lh_unknown']
    else: 
        atlas = 'HCPMMP1'
        _,_,rh_names = read_annot(f"../data/parcellations/rh.{atlas}.annot")
        drop = rh_names    
    
    spin_maps = nnsurf.spin_data(
        data = np.array(maps),
        drop = drop,
        version = "fsaverage",
        lhannot = f"../data/parcellations/lh.{atlas}.annot",
        rhannot = f"../data/parcellations/rh.{atlas}.annot",
        n_rotate = n
    )
    np.save(outfile, spin_maps)



from neuromaps.datasets import fetch_atlas
from neuromaps.nulls.spins import get_parcel_centroids, gen_spinsamples
def generate_spins(n=1000, blocks=1):
    """
    First generate spins. Do this once so it's not repeated for each gradient
    Generate in blocks of 1000 to prevent crashes
    """
    # Get sphere from which to generate spins (surface, not volume)
    surfaces = fetch_atlas('fsaverage', '41k')['sphere']
    # Create hemisphere ids needed to make spins
    # (This function is named 'get_parcel_centroids' but for cornblath spin method we use vertices)
    coords, hemiid = get_parcel_centroids(surfaces, method='surface')
    # Actually generate the spins
    for i in range(blocks):
        spins = gen_spinsamples(coords, hemiid, n_rotate=n, verbose=1)
        if blocks==0:
            np.save(f"../outputs/permutations/spins_41k_{n}.npy", spins)
        else:
            np.save(f"../outputs/permutations/spins_41k_{n}_{i}.npy", spins)
        print(f"\nGenerated block {i} of {n} spins")
    

from neuromaps.nulls import cornblath
def generate_nulls_from_gradients(scores, spins, hcp_img=None, n=10, n_components=3, only_left=True,
                           outfile='../outputs/permutations/spin_gradients_10.npy'):
    ## Next, get parcellation files in the same surface space as the spins (fsaverage)
    if hcp_img is None:
        hcp_img_files = ('../data/parcellations/lh.HCPMMP1.annot',
                         '../data/parcellations/rh.HCPMMP1.annot')
        hcp_img = annot_to_gifti(hcp_img_files)

    ## Reindex scores to have NA where parcels are missing (including all right hemi)
    scores_reindex = scores.reindex(range(1,361)).iloc[:,:n_components].values
    ### Drop parcels where data are missing in the 10k fsaverage HCPMMP parcellation template
    ## scores_reindex = np.delete(scores_reindex, [120,300], axis=0)

    ## Finally, for each gradient, compute nulls by projecting up to vertices and reaveraging
    null_scores = np.zeros([360, scores_reindex.shape[1], n])
    for i in range(3):
        _scores = scores_reindex[:,i]
        null_scores[:,i,:] = cornblath(
                data=_scores, 
                atlas='fsaverage', density='10k', 
                parcellation=hcp_img, 
                n_perm=n, spins=spins
                )
        print(f"\nGenerated {n} null models of axis {i}")
    
    # Drop right hemi
    if only_left:
        null_scores = null_scores[:180,:,:]

    np.save(outfile, null_scores)


# LEGACY VERSION
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

    
def generate_simulations(maps, n=10,
                         atlas = 'hcp',
                         dist_mat = None,
                         outfile='../outputs/sim_maps_1000.npy'):
    """
    Generate null maps using brainsmash
    """
    
    if 'label' in maps.columns:
        maps=maps.drop('label', axis=1)
    
    if atlas == 'hcp' and dist_mat is None:
        dist_mat="../data/parcellations/LeftParcelGeodesicDistmat.txt"
    elif atlas == 'dk':
        dist_mat="../data/parcellations/LeftParcelGeodesicDistmat_DK.txt"
    
    null_maps = np.zeros([maps.shape[0], maps.shape[1], n])
    
    for m in range(maps.shape[1]):
        base_map = maps.iloc[:,m].values
        base = Base(x=base_map, D=dist_mat)
        nulls = base(n)
        null_maps[:,m,:] = nulls.swapaxes(0,1)
        
    np.save(outfile, null_maps)


def corr_nulls_from_grads(null_grads, scores, maps, method='pearsonr', 
                          reindex=True, pool=False, pool_frac=.3, 
                          adjust='fdr_bh', adjust_by_label=False,
                          n_components=3):
    """
    Get correlations with maps from gradient score nulls
    Uses numpy masked array to handle missing values
    """
    # Filter maps for regions in gradients
    maps = maps.set_axis(range(1, maps.shape[0]+1)).loc[scores.index, :]
    # Filter nulls for regions in gradients
    # null_grads = null_grads[scores.index-1, :, :]

    # # Reindex scores to all regions
    # if reindex:
    #     scores = scores.reindex(range(1,181))
    # null_grads = null_grads

    n_maps = maps.shape[1]
    output_frame = np.zeros((n_components*n_maps, 2))
    output_index = pd.MultiIndex.from_product([scores.iloc[:,:n_components].columns, maps.columns])

    # For each gradient
    for g in range(null_grads.shape[1]):
        _scores = scores.iloc[:,g].values
        _scores = _scores.astype(np.longdouble) # force ValueError in compare_images
        # Optionally pool maps together
        if pool:
            # Set dimensions for reshaping
            # Downsample nulls because of pooling
            m = null_grads.shape[1]
            n = int(null_grads.shape[2] * pool_frac)
            null_grads_downsample = null_grads[:,:,:n]
            _nulls = null_grads_downsample.reshape(-1, g*n)
        else:
            _nulls = null_grads[:,g,:]
        
        # For each map
        for m in range(maps.shape[1]):
            _map = maps.iloc[:,m].values
            _map = _map.astype(np.longdouble) # force compare_images to work
            _r, _p = compare_images(_scores, _map, nulls=_nulls, metric=method, ignore_zero=False)
            output_frame[m+g*n_maps,:] = [_r, _p]

    # Output clean dataframe
    output_adjusted = (pd.DataFrame(output_frame, index=output_index) 
                        .set_axis(['r','p'], axis=1)
                        .rename_axis(['G','map']).reset_index()
    )
    
    # Multiple comparisons
    if adjust is not None:
        # Adjust across axes only (not by label)?
        if adjust_by_label:
            output_adjusted = (output_adjusted
                .assign(q = lambda x: x.groupby('G')
                                       .apply(lambda y: pd.Series(multipletests(y['p'], method=adjust)[1], index=y.index))
                                       .reset_index(0, drop=True) # make index match
                                       )
                )
        else:
            output_adjusted = (output_adjusted
                .assign(q = lambda x: multipletests(x['p'], method=adjust)[1])
                )
    else:       
        output_adjusted = output_adjusted.assign(q = lambda x: x['p'])

    return output_adjusted


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
        version_scores = version.clean_scores(n_components=n_components)
        corrs = get_corrs(version_scores, pca_scores, n_components=n_components).abs()
        matches = version.match_components(corrs, n_components=n_components)
        index_ = [f'MRI PC{i+1}' for i in range(n_components)]
        var = pd.Series(pca_var, index=index_, name='var')

        MRI_PC = [ix[0] for ix in matches[0]]
        G = [ix[1] for ix in matches[0]]
        pca_corrs_boot[i] = (pd.DataFrame({
                                'MRI_PC':MRI_PC, 
                                'G':G, 
                                'boot':i,
                                'r':matches[1]
                                })
                        .set_index('MRI_PC').join(var)
                        .reset_index()
                        )

    pca_corrs_boot = pd.concat(pca_corrs_boot).reset_index(drop=True)
    return pca_corrs_boot


# Atlas regions enrichment
def compute_region_means(scores, null_scores, regions_series, method='median'):
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
    
    true = pd.concat(true).unstack(1).set_axis(['G1','G2','G3'], axis=1)
    null = pd.concat(null).set_axis(['G1','G2','G3'], axis=1) \
                 .reset_index(level=0).rename({'level_0':'label'}, axis=1)

    return true, null




###
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