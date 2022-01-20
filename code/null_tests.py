# Code for comparing spin tests with surrogate map tests
import numpy as np, pandas as pd
from processing_helpers import *
from analysis_helpers import *

from netneurotools import freesurfer as nnsurf
from nibabel.freesurfer.io import read_annot
from brainsmash.mapgen.base import Base


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
    

def generate_spins_from_pcs(scores,
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
        n_rotate = 1000,
    )
    np.save(outfile, spin_pcs)

    
def generate_surrogates(maps, n=10,
                       outfile='../outputs/sim_maps_1000.npy'):
    """
    Generate null maps using brainsmash
    """
    dist_mat_file = "../data/LeftParcelGeodesicDistmat.txt"
    null_maps = np.zeros([maps.shape[0], maps.shape[1], n])
    
    for m in range(maps.shape[1]):
        base_map = maps.iloc[:,m].values
        base = Base(x=base_map, D=dist_mat_file)
        nulls = base(n)
        null_maps[:,m,:] = nulls.swapaxes(0,1)
        
    np.save(outfile, null_maps)
    
def corr_nulls_from_maps(null_maps, scores, maps, method='pearson'):
    """
    Get correlations with PC scores from null maps
    """
    null_corrs = {}
    for m, mapname in enumerate(maps.columns):
        nulls = pd.DataFrame(null_maps[:,m,:],
                             index=list(range(1,181)))
        null_corrs[mapname] = (
            pd.concat([scores, nulls], axis=1)
            .corr(method=method).iloc[3:,:3]
        )
        
    null_corrs = (
        pd.concat(null_corrs).reset_index(level=0)
        .set_axis(['map','PC1','PC2','PC3'],axis=1)
    )
    return null_corrs

def corr_nulls_from_pcs(null_pcs, scores, maps, method='pearson'):
    """
    Get correlations with maps from PC score nulls
    """
    n_maps = maps.shape[1]
    null_corrs = {}
    for m in range(null_pcs.shape[1]):
        nulls = pd.DataFrame(null_pcs[:,m,:],
                             index=scores.index)
        null_corrs[m] = (
            pd.concat([maps.set_axis(list(range(1,181))), nulls], axis=1)
            .corr(method=method).iloc[:n_maps,n_maps:]
            .stack().droplevel(1)
        )
        
    null_corrs = (
        pd.concat(null_corrs, axis=1).reset_index(level=0)
        .set_axis(['map','PC1','PC2','PC3'], axis=1)
    )
    return null_corrs



