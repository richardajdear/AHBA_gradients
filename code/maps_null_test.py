# Code for null tests (e.g. spin test) of spatial maps
import numpy as np, pandas as pd
from neuromaps.stats import compare_images
from neuromaps.datasets import fetch_atlas
from neuromaps.nulls.spins import get_parcel_centroids, gen_spinsamples
from neuromaps.nulls.spins import parcels_to_vertices, vertices_to_parcels
from neuromaps.nulls import cornblath
from neuromaps.transforms import fsaverage_to_fsaverage
from neuromaps.images import annot_to_gifti
from statsmodels.stats.multitest import multipletests
# from brainsmash.mapgen import Base


def correlate_maps_with_null_scores(
            null_scores, scores, maps, 
            method='pearsonr', reindex=True,
            adjust='fdr_bh', adjust_by_label=False
        ):
    """
    Spin test a set of maps against a set of scores with null scores
    Uses neuromaps.compare_images to perform the test
    """

    # Filter maps for regions in scores - assumes regions are in same order
    # maps = maps.set_axis(range(1, maps.shape[0]+1)).loc[scores.index, :]
    
    # Optional: filter nulls for regions in scores
    # null_scores = null_scores[scores.index-1, :, :]

    # Optional: reindex scores to all regions
    if reindex:
        scores = scores.reindex(range(1, null_scores.shape[0]+1))

    n_maps = maps.shape[1]
    n_components = scores.shape[1]
    output_frame = np.zeros((n_components*n_maps, 2))
    output_index = pd.MultiIndex.from_product([scores.iloc[:,:n_components].columns, maps.columns])

    # For each component:
    for C in range(null_scores.shape[1]):
        _scores = scores.iloc[:,C].values.astype(np.longdouble)
        _nulls = null_scores[:,C,:]
        
        # Test each map
        for m in range(maps.shape[1]):
            _map = maps.iloc[:,m].values.astype(np.longdouble)
            _r, _p = compare_images(_scores, _map, nulls=_nulls, metric=method, ignore_zero=False)
            output_frame[m+C*n_maps,:] = [_r, _p]

    # Output clean dataframe
    output = (pd.DataFrame(output_frame, index=output_index) 
                        .set_axis(['r','p'], axis=1)
                        .rename_axis(['C','map']).reset_index()
    )
    
    # Multiple comparisons
    if adjust is not None:
        # Adjust for only across the comparisons within each component
        if adjust_by_label:
            output_adjusted = (output
                .assign(q = lambda x: x.groupby('C')
                                       .apply(lambda y: pd.Series(multipletests(y['p'], method=adjust)[1], index=y.index))
                                       .reset_index(0, drop=True) # make index match
                                       )
                )
        # Or adjust across all comparisons for all components (more strict)
        else:
            output_adjusted = (output
                .assign(q = lambda x: multipletests(x['p'], method=adjust)[1])
                )
    else:       
        # Or don't adjust at all
        output_adjusted = output.assign(q = lambda x: x['p'])

    return output_adjusted



def generate_shuffles(maps, n=1000, outfile='../outputs/shuffle_maps_1000.npy'):
    """
    Random shuffles null model for comparison
    """
    if 'label' in maps.columns:
        maps = maps.drop('label', axis=1)
    
    shuffle_maps = np.repeat(maps.values[:,:,np.newaxis], 1000, axis=2)
    for i in range(shuffle_maps.shape[2]):
        np.random.shuffle(shuffle_maps[:,:,i])

    np.save(outfile, shuffle_maps)



def generate_spins(n=1000, blocks=1, density='41k', 
                   save_name=None, save_dir="../outputs/permutations"):
    """
    Generate spins of the fsaverage sphere
    (This saves time as we can use the same spins to make nulls for multiple maps)
    Optionally generate in multiple blocks (helps if getting out of memory errors)
    """
    # Get sphere from which to generate spins (surface, not volume)
    surfaces = fetch_atlas('fsaverage', density)['sphere']
    # Create hemisphere ids needed to make spins
    # (This function is named 'get_parcel_centroids' but for cornblath spin method we use vertices)
    coords, hemiid = get_parcel_centroids(surfaces, method='surface')

    # Auto-generate a name for the saved spins if not given
    if save_name is None:
        save_name = f"spins_{density}_{n}"

    # Generate the spins
    for i in range(blocks):
        spins = gen_spinsamples(coords, hemiid, n_rotate=n, verbose=1)
        if blocks==1:
            save_path = f"{save_dir}/{save_name}.npy"
            np.save(save_path, spins)
            print(f"\nSaved spins to {save_path}")
        else:
            save_path = f"{save_dir}/{save_name}_{i}.npy"
            np.save(save_path, spins)
            print(f"\nSaved block {i} of spins to {save_path}")
    print(f"\nGenerated {blocks} blocks of {n} spins at density {density}")
    


def generate_nulls_from_components(scores, spins, atlas='hcp', parcellation_img=None, density='41k', 
                                  n=10, only_left=True,
                                  save_dir = '../outputs/permutations/',
                                  save_name = 'spin_10'):
    """
    Generate null models using predefined spins
    Use fsaverage 41k (10k has an error)
    """
    ## First, get parcellation files in the same surface space as the spins (fsaverage)
    if parcellation_img is None and atlas=='hcp':
        parcellation_img_files = ('../data/parcellations/lh.HCPMMP1.annot',
                                  '../data/parcellations/rh.HCPMMP1.annot')
        parcellation_img = annot_to_gifti(parcellation_img_files)
        parcellation_img = fsaverage_to_fsaverage(parcellation_img, target_density=density, method='nearest')
    elif parcellation_img is None and atlas=='dk':
        parcellation_img_files = ("../data/parcellations/lh.aparc.annot", 
                                  "../data/parcellations/rh.aparc.annot")
        parcellation_img = annot_to_gifti(parcellation_img_files)
        parcellation_img = fsaverage_to_fsaverage(parcellation_img, target_density=density, method='nearest')

    n_components = scores.shape[1]

    ## Reindex scores to have NA where parcels are missing (including all right hemi)
    if atlas=='hcp':
        scores_reindex = scores.reindex(range(1,361)).iloc[:,:n_components].values
        # if density=='10k':
        # ### Drop parcels where data are missing in the 10k fsaverage HCPMMP parcellation template
        #     scores_reindex = np.delete(scores_reindex, [120,300], axis=0)
    elif atlas=='dk':
        scores_reindex = scores.reindex(range(1,69)).iloc[:,:n_components].values


    ## Finally, for each component, compute nulls by projecting up to vertices and reaveraging
    ## Uses cornblath method; can replace with other methods from https://netneurolab.github.io/neuromaps/api.html#module-neuromaps.nulls
    null_scores = np.zeros([scores_reindex.shape[0], n_components, n])
    for i in range(n_components):
        _scores = scores_reindex[:,i]
        null_scores[:,i,:] = cornblath(
                data = _scores,
                atlas = 'fsaverage', 
                density = density, 
                parcellation = parcellation_img, 
                n_perm = n, 
                spins = spins
            )
        print(f"\nGenerated {n} null spins of component {i}")
    
    # Drop right hemi
    if only_left and atlas=='hcp':
        null_scores = null_scores[:180,:,:]
    elif only_left and atlas=='dk':
        null_scores = null_scores[:34,:,:]

    save_path = f"{save_dir}{save_name}.npy"
    np.save(save_path, null_scores)
    print(f"Saved null spins to {save_path}")


def generate_simulations(maps, n=10,
                         atlas = 'hcp',
                         dist_mat = None,
                         outfile='../outputs/sim_maps_1000.npy'):
    """
    Alternate null method: Generate null maps using brainsmash
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




###

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
