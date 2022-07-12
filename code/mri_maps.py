# Code for comparing PC maps with MRI maps
import numpy as np, pandas as pd
from processing_helpers import *
from analysis_helpers import *

from netneurotools import freesurfer as nnsurf
from nibabel.freesurfer.io import read_annot
from brainsmash.mapgen.base import Base
from scipy.stats import percentileofscore
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests
from nilearn.input_data import NiftiLabelsMasker




def get_maps(data_dir="../data/stat_maps_HCP_forRichard.csv", filter=True):
    maps = (
        pd.read_csv(data_dir)
        .apply(lambda x: (x-np.mean(x))/np.std(x))
        .sort_index(axis=1)
        .set_index(get_labels_hcp()[:180])
        # .rename_axis('region')#.reset_index()
    )
    
    selected_maps = [
        'T1T2',
        'thickness',
        #'glasser_CMRO2',
        #'glasser_CMRGlu',        
        'G1_fMRI',        
        'PC1_neurosynth',        
        'externopyramidisation',        
        'glasser_GI',
        'hill.evo_remapped',        
        'hill.dev_remapped',
        'glasser_CBF',
        'allom',
        
        # 'glasser_CBV',
        # 'asl',
        
        # 'x':'X axis',
        # 'y':'Y axis',
        # 'z':'Z axis'
    ]
    
    if filter:
        maps = maps.loc[:, selected_maps]
        
    return maps


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


def get_corrs(scores, maps, method='pearson', atlas='hcp'):
    
    # labels = get_labels_dk() if atlas=='dk' else get_labels_hcp()[:180]
        
    corrs = (
        scores
        .set_index('label')
        .join(maps)
        .corr(method=method).iloc[3:,:3]
        .set_axis([f'G{i+1}' for i in range(3)], axis=1)
    )
    return corrs


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

def generate_spins_from_gradients(scores, n=10,
                           outfile='../outputs/permutations/spin_gradients_1000.npy',
                           atlas='hcp'):
    
    if atlas == 'dk':
        atlas = 'aparc'
        _,_,rh_names = read_annot(f"../data/parcellations/rh.{atlas}.annot")
        drop = rh_names + ['lh_corpuscallosum', 'lh_unknown']
    else: 
        atlas = 'HCPMMP1'
        _,_,rh_names = read_annot("../data/parcellations/rh.HCPMMP1.annot")
        # Find regions missing in scores
        missing_rois = list(set(get_labels_hcp()[:180]).difference(scores['label']))
        # Fix 7Pl to match HCP codes
        missing_rois = ['7PL' if x=='7Pl' else x for x in missing_rois]
        # Drop missing regions
        rh_names = rh_names + [np.bytes_('L_' + roi + '_ROI') for roi in missing_rois]
        drop = rh_names
    
    spin_pcs = nnsurf.spin_data(
        data = np.array(scores.set_index('label')),
        drop = drop,
        version = "fsaverage",
        lhannot = f"../data/parcellations/lh.{atlas}.annot",
        rhannot = f"../data/parcellations/rh.{atlas}.annot",
        n_rotate = n,
    )
    np.save(outfile, spin_pcs)

    
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
        dist_mat="../data/LeftParcelGeodesicDistmat.txt"
    elif atlas == 'dk':
        dist_mat="../data/LeftParcelGeodesicDistmat_DK.txt"
    
    null_maps = np.zeros([maps.shape[0], maps.shape[1], n])
    
    for m in range(maps.shape[1]):
        base_map = maps.iloc[:,m].values
        base = Base(x=base_map, D=dist_mat)
        nulls = base(n)
        null_maps[:,m,:] = nulls.swapaxes(0,1)
        
    np.save(outfile, null_maps)

    
def corr_nulls_from_maps(null_maps, scores, maps, method='pearson', pool=False, pool_frac=.1):
    """
    Get correlations with scores from null maps
    """
    
    # Drpp 'label' field if it exists
    if 'label' in scores.columns:
        scores = scores.drop('label', axis=1)
        
    if pool:
        # Set dimensions for reshaping
        # Downsample nulls because of pooling
        m = null_maps.shape[1]
        n = int(null_maps.shape[2] * pool_frac)
        null_maps_downsample = null_maps[:,:,:n]
        null_maps_pool = null_maps_downsample.reshape(-1, m*n)
    
    null_corrs = {}
    # Define index as 1 to n
    index = [i+1 for i in range(null_maps.shape[0])]
    for m, mapname in enumerate(maps.columns):
        # Optionally pool maps together
        if pool:
            nulls = pd.DataFrame(null_maps_pool, index=index)
        else:
            nulls = pd.DataFrame(null_maps[:,m,:], index=index)
        
        # Concat scores and nulls, matching on index directly
        df_concat = pd.concat([scores, nulls], axis=1)
        # Optionally correlate with spearman, otherwise pearson
        if method == 'spearman':
            # .rank().corr() is faster than .corr(method='spearman') because of null checks
            df_corr = df_concat.rank().corr()
        else:
            df_corr = df_concat.corr()
        # Cleanup
        null_corrs[mapname] = df_corr.iloc[3:,:3]
        
    # Concat rows (each map)
    null_corrs = (
        pd.concat(null_corrs).reset_index(level=0)
        .set_axis(['map','G1','G2','G3'],axis=1)
    )
    return null_corrs


def corr_nulls_from_grads(null_grads, scores, maps, method='pearson', pool=False, pool_frac=.1):
    """
    Get correlations with maps from gradient score nulls
    """
    if pool:
        # Set dimensions for reshaping
        # Downsample nulls because of pooling
        m = null_grads.shape[1]
        n = int(null_grads.shape[2] * pool_frac)
        null_grads_downsample = null_grads[:,:,:n]
        null_grads_pool = null_grads_downsample.reshape(-1, m*n)
    
    null_corrs = {}
    # Gradients do not have all regions, so take scores index as index
    index = scores.index
    # Maps still don't have an index, so add it here
    maps_index = [i+1 for i in range(maps.shape[0])]
    for m in range(null_grads.shape[1]):
        # Optionally pool maps together
        if pool:
            nulls = pd.DataFrame(null_grads_pool, index=index)
        else:
            nulls = pd.DataFrame(null_grads[:,m,:], index=index)
        
        # Concat scores and nulls
        df_concat = pd.concat([maps.set_axis(maps_index), nulls], axis=1)
        # Optionally correlate with spearman, otherwise pearson
        if method == 'spearman':
            # .rank().corr() is faster than .corr(method='spearman') because of null checks
            df_corr = df_concat.rank().corr()
        else:
            df_corr = df_concat.corr()
        # Cleanup and stack maps into vector
        n_maps = maps.shape[1]
        null_corrs[m] = df_corr.iloc[:n_maps,n_maps:].stack().droplevel(1)
        
    # Concat columns (each gradient)
    null_corrs = (
        pd.concat(null_corrs, axis=1).reset_index(level=0)
        .set_axis(['map','G1','G2','G3'], axis=1)
    )
    return null_corrs
    
    


def get_null_p(true_corrs, null_corrs, adjust='fdr_bh', order=None):
    """
    Get p values
    """
    # For each map, for each set of nulls, compute the percentile of the true correlations
    null_pct = np.zeros(true_corrs.shape)
    for m, _map in enumerate(true_corrs.index):
        for i in range(3):
            _null_corrs = null_corrs.set_index('map').loc[_map].iloc[:,i]
            _corr = true_corrs.iloc[m,i]
            pct = percentileofscore(_null_corrs, _corr)/100
            null_pct[m,i] = pct

    true_mean = true_corrs.stack().rename('true_mean')

    null_mean = (null_corrs
                 .groupby('map').agg(['mean','std'])
                 .stack(0)
                 .rename_axis([None,None])
                 .set_axis(['null_mean', 'null_std'], axis=1)
                )
            
    null_p = (pd.DataFrame(null_pct, index=true_corrs.index, columns=true_corrs.columns)
              .stack().rename('pct').to_frame()
              .join(true_mean)
              .join(null_mean)
              .assign(z = lambda x: (x['true_mean'] - x['null_mean'])/x['null_std'])
              .assign(pos = lambda x: [pct > 0.5 for pct in x['pct']])
              .assign(p = lambda x: [(1-pct)*2 if pct>0.5 else pct*2 for pct in x['pct']]) # x2 for bidirectional
              .assign(sig = lambda x: x['p'] < .05)
             )
    
    if adjust is not None:
        null_p = (null_p
             .assign(q = lambda x: multipletests(x['p'], method=adjust)[1])
             .assign(sig = lambda x: x['q'] < .05)
             # .assign(q_abs = lambda x: [1-q if pos else q for pos, q in zip(x['pos'], x['q'])])
            )
    else:
        null_p = (null_p
             .assign(q = lambda x: x['p'])
             .assign(sig = lambda x: x['p'] < .05)
                 )
    
    null_p = (null_p
              .reset_index()
              .rename({'level_0':'map', 'level_1':'G'}, axis=1)
              # .assign(map = lambda x: pd.Categorical(x['map'], ordered=True, categories=order))
             )
    
    # Fix order of maps
    if order is None:
        order = true_corrs.index
    null_p = (null_p
              .assign(map = lambda x: pd.Categorical(x['map'], ordered=True, categories=order))
              .sort_values('map')
         )

    return null_p



def maps_pca(maps, short_names=None):
    """
    Run PCA on maps
    """
    cols = ['MRI PC1','MRI PC2','MRI PC3']
    pca = PCA(n_components=3).fit(maps)
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
    
    















    

def print_corrs_sig(corrs, null_p):
    map_corrs_sig = (corrs
     .loc[null_p.index]
     .round(2).astype('string')
     .where(null_p > .05, other = lambda x: x+' *')
     .where(null_p > .01, other = lambda x: x+'*')
     .where(null_p > .001, other = lambda x: x+'*')
    )
    # map_corrs_sig.to_csv("../outputs/map_corrs_sig.csv")
    return map_corrs_sig


def get_inflammation_data(img="../data/inflammation/Activation_proportion.nii.gz"):
    hcp_mask = "../data/parcellations/HCP-MMP_1mm.nii.gz"
    mask = NiftiLabelsMasker(hcp_mask, resampling_target='data')
    inflammation_img = img
    inflammation_hcp = mask.fit_transform(inflammation_img).squeeze()

    inflammation_hcp = (
        pd.Series(inflammation_hcp)
        .rename('activation')
        .to_frame()
        .set_axis([i+1 for i in range(360)])
        .join(get_labels_hcp())
        .set_index('label')
        # .rename({'label':'region'}, axis=1)
        # set_index('region')
        # .assign(hemi=lambda x:['left' if id<=180 else 'right' for id in x.index])
               )
    return inflammation_hcp