# Code for comparing PC maps with MRI maps
import numpy as np, pandas as pd
from processing_helpers import *
from analysis_helpers import *

from netneurotools import freesurfer as nnsurf
from nibabel.freesurfer.io import read_annot
from brainsmash.mapgen.base import Base
from scipy.stats import percentileofscore




def get_maps(data_dir="../data/stat_maps_HCP_forRichard.csv"):
    maps_jakob = (
        pd.read_csv("../data/stat_maps_HCP_forRichard.csv")
        .apply(lambda x: (x-np.mean(x))/np.std(x))
        .sort_index(axis=1)
        .set_index(get_labels_hcp()[:180])
        .rename_axis('region')#.reset_index()
    )
    
    selected_maps = {
        'T1T2':'Myelination (T1w/T2w ratio)',
        'thickness':'Cortical Thickness',    
        'externopyramidisation':'Intracortical Connectivity Distance',
        'allom':'Size Variation (allometric scaling)',
        # 'x':'X axis',

        'CBF':'Cerebral Blood Flow',
        'glasser_CMRGlu':'Glucose Metabolic Rate (CMRGlu)',
        'glasser_CMRO2':'Oxygen Metabolic Rate (CMRO2)',
        'glasser_GI':'Glycolytic Index (CMRO2/CMRGlu)',
        # 'y':'Y axis',

        'G1_fMRI': 'fMRI PC1',
        'PC1_neurosynth': 'NeuroSynth PC1',
        'hill.evo_remapped':'Evolutionary Expansion',
        'hill.dev_remapped':'Developmental Expansion',
        # 'z':'Z axis'
    }

    maps = maps_jakob.loc[:, list(selected_maps.keys())].set_axis(selected_maps.values(), axis=1)
    return maps


def get_corrs(scores, maps, method='pearson'):
    corrs = (
        scores
        .join(get_labels_hcp()[:180])
        .set_index('label')
        .join(maps) #.set_index('region'))
        .corr(method=method).iloc[3:,:3]
        .set_axis([f'PC{i+1}' for i in range(3)], axis=1)
    )
    return corrs


def generate_shuffles(maps, n_shuffles=1000,
                      outfile='../outputs/shuffle_maps_1000.npy'):
    shuffle_maps = np.repeat(maps.values[:,:,np.newaxis], 1000, axis=2)
    for i in range(shuffle_maps.shape[2]):
        np.random.shuffle(shuffle_maps[:,:,i])

    np.save(outfile, shuffle_maps)


def generate_spins(maps, n_rotate=1000, 
                   outfile='../outputs/spin_maps_1000.npy'):
    """
    Generate spins from a set of maps and save to file
    """
    _,_,rh_names = read_annot("../data/rh.HCPMMP1.annot")

    spin_maps = nnsurf.spin_data(
        data = np.array(maps),
        drop = rh_names,
        version = "fsaverage",
        lhannot = "../data/lh.HCPMMP1.annot",
        rhannot = "../data/rh.HCPMMP1.annot",
        n_rotate = n_rotate
    )
    np.save(outfile, spin_maps)

def generate_spins_from_pcs(scores, n_rotate=1000,
    outfile='../outputs/spin_pcs_1000.npy'):
    
    scores = scores.join(get_labels_hcp())
    
    _,_,rh_names = read_annot("../data/rh.HCPMMP1.annot")
    # Drop regions missing in PCs
    missing_rois = list(set(get_labels_hcp()[:180]).difference(scores['label']))
    rh_names = rh_names + [np.bytes_('L_' + roi + '_ROI') for roi in missing_rois]

    spin_pcs = nnsurf.spin_data(
        data = np.array(scores.set_index('label')),
        drop = rh_names,
        version = "fsaverage",
        lhannot = "../data/lh.HCPMMP1.annot",
        rhannot = "../data/rh.HCPMMP1.annot",
        n_rotate = n_rotate,
    )
    np.save(outfile, spin_pcs)

    
def generate_surrogates(maps, n=10,
                        dist_mat="../data/LeftParcelGeodesicDistmat.txt",
                        outfile='../outputs/sim_maps_1000.npy'):
    """
    Generate null maps using brainsmash
    """
    null_maps = np.zeros([maps.shape[0], maps.shape[1], n])
    
    for m in range(maps.shape[1]):
        base_map = maps.iloc[:,m].values
        base = Base(x=base_map, D=dist_mat)
        nulls = base(n)
        null_maps[:,m,:] = nulls.swapaxes(0,1)
        
    np.save(outfile, null_maps)

    
def corr_nulls_from_maps(null_maps, scores, maps, method='pearson', pool=False):
    """
    Get correlations with PC scores from null maps x
    """
    if pool:
        mxn = null_maps.shape[1] * null_maps.shape[2] # 2nd dimension of array after pooling, for rehsaping
        null_maps_pool = null_maps.reshape(-1, mxn)
    
    null_corrs = {}
    for m, mapname in enumerate(maps.columns):
        # Optionally pool maps together
        if pool:
            nulls = pd.DataFrame(null_maps_pool, index=list(range(1,181)))
        else:
            nulls = pd.DataFrame(null_maps[:,m,:], index=list(range(1,181)))
        
        # Concat PC scores and nulls
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
        .set_axis(['map','PC1','PC2','PC3'],axis=1)
    )
    return null_corrs


def corr_nulls_from_pcs(null_pcs, scores, maps, method='pearson', pool=False):
    """
    Get correlations with maps from PC score nulls
    """
    if pool:
        mxn = null_pcs.shape[1] * null_pcs.shape[2] # 2nd dimension of array after pooling, for rehsaping
        null_pcs_pool = null_pcs.reshape(-1, mxn)
    
    null_corrs = {}
    for m in range(null_pcs.shape[1]):
        # Optionally pool maps together
        if pool:
            nulls = pd.DataFrame(null_pcs_pool, index=scores.index)
        else:
            nulls = pd.DataFrame(null_pcs[:,m,:], index=scores.index)
        
        # Concat PC scores and nulls
        df_concat = pd.concat([maps.set_axis(list(range(1,181))), nulls], axis=1)
        # Optionally correlate with spearman, otherwise pearson
        if method == 'spearman':
            # .rank().corr() is faster than .corr(method='spearman') because of null checks
            df_corr = df_concat.rank().corr()
        else:
            df_corr = df_concat.corr()
        # Cleanup and stack maps into vector
        n_maps = maps.shape[1]
        null_corrs[m] = df_corr.iloc[:n_maps,n_maps:].stack().droplevel(1)
        
    # Concat columns (each PC)
    null_corrs = (
        pd.concat(null_corrs, axis=1).reset_index(level=0)
        .set_axis(['map','PC1','PC2','PC3'], axis=1)
    )
    return null_corrs

    'spin_maps_p', 'spin_maps_p_pool',
    'sim_maps_p', 'sim_maps_p_pool',
    'spin_pcs_p', 'spin_pcs_p_pool',
    'sim_pcs_p', 'sim_pcs_p_pool',
    
    


def get_null_p(corrs, null_corrs):
    """
    Get p values
    """
    null_p = np.zeros(corrs.shape)
    for m, _map in enumerate(corrs.index):
        for i in range(3):
            _null_corrs = null_corrs.set_index('map').loc[_map].iloc[:,i]
            _corr = corrs.iloc[m,i]
            p = percentileofscore(_null_corrs, _corr)/100
            if p > .5:
                p = 1-p
            null_p[m,i] = p

    null_p = pd.DataFrame(null_p, index=corrs.index, columns=corrs.columns)
    
    return null_p
    
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