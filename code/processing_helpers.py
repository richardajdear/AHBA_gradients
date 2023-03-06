# Helper functions for processing

import numpy as np, pandas as pd
import pickle
import abagen
import abagen_allen_tweaked
from abagen.correct import keep_stable_genes
from abagen.utils import efficient_corr
from itertools import combinations
import warnings
import nibabel as nib
from neuromaps.images import annot_to_gifti
from neuromaps.nulls.spins import parcels_to_vertices, vertices_to_parcels




def save_pickle(data, fname):
    with open('../outputs/' + fname + '.pickle', 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_pickle(fname):
    with open('../outputs/' + fname + '.pickle', 'rb') as handle:
        return pickle.load(handle)



def get_expression_abagen(atlas, 
                          save_name=None,
                          data_dir='../data/abagen-data',
                          return_donors=False,
                          return_counts=False,
                          return_labels=False,
                          return_stability=False,
                          verbose=0,
                          only_left=True,
                          only_cortex=True,
                          ibf_threshold=0.5,
                          probe_selection='diff_stability',
                          lr_mirror='rightleft',
                          region_agg='donors',
                          tolerance=2,
                          sample_norm='srs',
                          gene_norm='srs',
                          donors='all',
                          donors_threshold=0,
                          gene_stability_threshold=0,
                          region_stability_threshold=0,                          
                          **kwargs):
    """
    Get expression matrix and matching labels from Abagen
    """
    
    # Pack overridden defaults into kwargs to pass to abagen
    kwargs = dict(verbose=verbose,
                  only_left=only_left,
                  only_cortex=only_cortex,
                  ibf_threshold=ibf_threshold,
                  probe_selection=probe_selection,
                  lr_mirror=lr_mirror,
                  region_agg=region_agg,
                  tolerance=tolerance,
                  sample_norm=sample_norm,
                  gene_norm=gene_norm,
                  donors=donors,
                  donors_threshold=donors_threshold,
                  **kwargs)

    with warnings.catch_warnings(): 
        warnings.simplefilter('ignore')
        _out = abagen_allen_tweaked.get_expression_data(
            atlas=atlas['image'],
            atlas_info=atlas['info'],
            data_dir=data_dir + '/microarray',
            n_proc=32,
            return_counts=return_counts,
            return_labels=return_labels,
            return_donors=True, # always true to get DS
            **kwargs)
        
    # If returning labels or counts, expression is the first element of tuple
    if return_labels or return_counts:
        expression_all = _out[0]
    else:
        expression_all = _out

    # Drop rows with all 0
    # expression_all = [e for e in expression_all]
    
    # Filter for stable regions
    # if region_stability_threshold > 0:
    #     expression_stable_regions, region_stability = keep_stable_regions(
    #         expression_all,
    #         threshold=region_stability_threshold,
    #         return_stability=True)
    #     region_stability = pd.Series(region_stability, index=expression_all[0].index)
    #     print(f'{expression_stable_regions[1].shape[1]} regions remain \
    #           after filtering for top {round(1-region_stability_threshold,2)} stability')

    # Filter for stable genes
    if gene_stability_threshold > 0:
        expression_stable_genes, gene_stability = keep_stable_genes(
            expression_all, 
            threshold=gene_stability_threshold,
            return_stability=True)
        gene_stability = pd.Series(gene_stability, index=expression_all[0].columns)
        print(f'{expression_stable_genes[0].shape[1]} genes remain \
              after filtering for top {round(1-gene_stability_threshold,2)} stability')

    # Combine donors together after filtering
    if return_donors:
        expression = expression_stable_genes
    else:
        expression = pd.concat(expression_stable_genes).groupby(level=0).mean()
    
    # Save combined expression data
    if save_name is not None:
        expression.to_csv(data_dir + "/expression/" + save_name + ".csv")
        if return_labels:
            _out[1].to_csv(data_dir + "/labels/" + save_name + ".csv")
    
    # Pack outputs into tuple
    out = (expression,)
    if return_labels:
        out += (_out[1],)
    if return_counts:
        if return_labels:
            out += (_out[2],)
        else:
            out += (_out[1],)
    if return_stability:
        out += (gene_stability, region_stability,)
    if len(out) == 1:
        out = out[0]
    
    return out


def keep_stable_regions(expression, threshold, 
                        percentile=True, rank=True, 
                        return_stability=True):
    """
    Keep only stable regions, using same logic as gene
    differential stability filter.
    """

    num_subj = len(expression)
    num_region = expression[0].shape[0]

    # flip axes to get stability of regions, not genes
    for_corr = expression.T if not rank else [e.T.rank() for e in expression]

    # for all donor pairs
    region_corrs = np.zeros((num_region, sum(range(num_subj))))
    for n, (s1, s2) in enumerate(combinations(range(num_subj), 2)):
        genes = np.intersect1d(for_corr[s1].dropna(axis=0, how='all').index,
                               for_corr[s2].dropna(axis=0, how='all').index)
        region_corrs[:, n] = efficient_corr(for_corr[s1].loc[genes],
                                            for_corr[s2].loc[genes])

    # average similarity across donors (ignore NaNs)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning,
                                message='Mean of empty slice')
        region_corrs = np.nanmean(region_corrs, axis=1)

    # calculate absolute threshold if percentile is desired
    if percentile:
        threshold = np.nanpercentile(region_corrs, threshold * 100)
    keep_regions = region_corrs >= threshold
    expression = [e.iloc[keep_regions, :] for e in expression]

    if return_stability:
        return expression, region_corrs

    return expression


    

def fetch_dk(native=True, only_cortex=True, only_left=False):
    """
    Get DK atlas and remove subcortex
    """
    atlas_dk = abagen.fetch_desikan_killiany(native=native)
    
    def remove_subcortex(img):
        if only_left:
            img_cortex = np.where((img.dataobj[:] > 34), 0, img.dataobj[:])
        else:
            img_cortex = np.where((img.dataobj[:] > 34) & (img.dataobj[:] < 48), 0, img.dataobj[:])
        return img.__class__(img_cortex, img.affine, img.header)
        
    if only_cortex:
        if native:
            for donor, img in atlas_dk['image'].items():
                img = nib.load(img)
                atlas_dk['image'][donor] = remove_subcortex(img)
    #         atlas_dk['info'] = pd.read_csv(atlas_dk['info']).replace('subcortex/brainstem','other')
        else:
            img = nib.load(atlas_dk['image'])
            atlas_dk['image'] = remove_subcortex(img)
    
    atlas_dk['info'] = pd.read_csv(atlas_dk['info'])
    
    return atlas_dk



def fetch_hcp(native=True, only_left=False):
    """
    Get HCP atlas
    """
    hcp_info = (
        pd.read_csv('../data/parcellations/HCP-MMP1_UniqueRegionList.txt')
        .rename(columns={'LR':'hemisphere', 'regionID':'id', 'region':'label'})
        .assign(structure='cortex')
        .assign(id=lambda x: [id-20 if id>180 else id for id in x['id']])
    )

    path = '../data/AHBAprocessing-data/parcellations/'
    name = 'HCPMMP1_acpc_uncorr.nii'
    donors = ['9861', '10021', '12876', '14380', '15496', '15697']
    
    if native:
        atlas_hcp = {
            'image':{donor:path + 'S%s/' % str(i+1) + name for i,donor in enumerate(donors)},
            'info':hcp_info
        }
    else:
        atlas_hcp = {
            'image':nib.load('../data/parcellations/HCP-MMP_1mm.nii.gz'),
            'info':pd.read_csv('../data/parcellations/HCP-MMP1_UniqueRegionList.txt').rename(columns={'LR':'hemisphere', 'regionID':'id', 'region':'label'}).assign(structure='cortex')
        }

    def remove_right(img):
        img_cortex = np.where((img.dataobj[:] > 180), 0, img.dataobj[:])
        return img.__class__(img_cortex, img.affine, img.header)
    
    if only_left:
        for donor, image in atlas_hcp['image'].items():
            atlas_hcp['image'][donor] = remove_right(image)

    return atlas_hcp


def fetch_hcp_ex():
    """
    Get HCP atlas with subcortex, data from https://github.com/wayalan/HCPex
    """

    atlas_hcp_ex = {
        'image': '../data/parcellations/HCPex.nii.gz',
        'info': '../data/parcellations/HCPex_info.csv'
    }

    return atlas_hcp_ex



def get_labels_dk():
    """
    Get DK atlas labels using Abagen and format for ggseg
    """
    atlas_dk = abagen.fetch_desikan_killiany(native=True)
    labels_dk = (
        pd.read_csv(atlas_dk['info'])
        .set_index('id')
        .assign(label = lambda x: 'lh_' + x['label'])['label']
    )
    return labels_dk


def get_labels_hcp():
    """
    Get HCP atlas labels from source file and format for ggseg
    """
    hcp_info = (
        pd.read_csv('../data/parcellations/HCP-MMP1_UniqueRegionList.txt')
        .rename(columns={'LR':'hemisphere', 'regionID':'id', 'region':'label'})
        .assign(structure='cortex')
        .assign(id=lambda x: [id-20 if id>180 else id for id in x['id']]) 
        # Aurina's HCP images code regions as 1-360, not 1-180,201-380
        # so recode to match the image files
    )
    labels_hcp = hcp_info.set_index('id')['label']
    return labels_hcp
    


def project_to_dk(scores,
                hcp_img_path = "../data/parcellations/lh.HCPMMP1.annot",
                dk_img_path = "../data/parcellations/lh.aparc.annot"
               ):
    """
    Project HCP scores to DK using annot files
    """
    hcp_img = annot_to_gifti(hcp_img_path)
    dk_img = annot_to_gifti(dk_img_path)

    scores_dk = np.zeros((34,3))
    for i in range(3):
        # Re-index gradient null values with NA
        g_hcp = scores[i].reindex(range(1,181)).values
        # Use HCP parcellation image to project HCP data to fsaverage
        g_fsaverage = parcels_to_vertices(g_hcp, hcp_img)
        # Use DK parcellation image to project fsaverage data into DK
        g_dk = vertices_to_parcels(g_fsaverage, dk_img)
        # Add to outputs
        scores_dk[:,i] = g_dk

    # Convert to dataframe
    scores_dk = pd.DataFrame.from_records(scores_dk, index=list(range(1,35)))

    return scores_dk

    

def fetch_dx():
    """
    Get Desterieux atlas
    """
    dx_info = (
        pd.read_csv("../data/parcellations/Destrieux.csv")
        .set_axis(['id','label'],axis=1)
        .assign(structure='cortex', hemisphere='L')
    )

    atlas_dx = {
        'image':'../data/parcellations/Destrieux_space-MNI152NLin6_res-4x4x4.nii.gz',
        'info': dx_info
    }

    return atlas_dx

def get_labels_dx():
    labels_dx = (
        pd.read_csv("../data/parcellations/desterieux_ggseg_labels.csv")
        .set_axis(['id','label'],axis=1)
        .set_index('id')['label']
    )
    return labels_dx


def fetch_schaefer(size=400):
    """
    Get Schaefer atlas
    """
    # s200_info = (
    #     pd.read_csv("../data/parcellations/schaefer200_ggseg_labels.csv")
    #     .assign(structure='cortex', hemisphere='L')
    #     #.loc[lambda x: x['id']!=1] # drop medial wall
    # )

    # Path to Schaefer200 image
    img_path = f'../data/parcellations/Schaefer2018_{size}Parcels_17Networks_order_FSLMNI152_1mm.nii.gz'
    
    # Drop right hemi
    img = nib.load(img_path)
    img_left = np.where((img.dataobj[:] > size/2+1), 0, img.dataobj[:])

    img_left = np.where((img_left == 1), 0, img_left)    # drop medial wall

    img = img.__class__(img_left, img.affine, img.header)
    
    atlas = {
        'image':img,
        'info': None
        #'info': s200_info
    }

    return atlas

def get_labels_schaefer(size=400):
    labels = (
        pd.read_csv(f"../data/parcellations/schaefer{size}_ggseg_labels.csv")
        .set_index('id')['label']
        .astype('str')
    )
    return labels

    


### Legacy



def get_labels_aseg():
    """
    Get aseg atlas labels using Abagen and format for ggseg
    """
    atlas_dk = abagen.fetch_desikan_killiany(native=True)
    labels_aseg = (
        pd.read_csv(atlas_dk['info'])
#         .query("structure != 'cortex'")
        .set_index('id')
        .assign(hemi_label = lambda x: x['hemisphere'].replace({'L':'Left-', 'R':'Right-', 'B':''}))
        .assign(label = lambda x: (x['label']
                                   .str.capitalize()
                                   .replace({'Thalamusproper':'Thalamus-Proper'}))
               )
        .assign(label = lambda x: x['hemi_label'] + x['label'])['label']
        .replace({'Brainstem':'brain-stem',
                  'Left-Rostralanteriorcingulate':'cc-anterior',
                  'Right-Rostralanteriorcingulate':'cc-anterior',
                  'Left-Caudalanteriorcingulate':'cc-mid-anterior',
                  'Right-Caudalanteriorcingulate':'cc-mid-anterior',
                  'Left-Posteriorcingulate':'cc-mid-posterior',
                  'Right-Posteriorcingulate':'cc-mid-posterior',
                  'Left-Isthmuscingulate':'cc-posterior',
                  'Right-Isthmuscingulate':'cc-posterior',
                 })
    )
    return labels_aseg



def read_expression_aurina(path):
    """
    Get expression matrix from Aurina
    """
    genes = pd.read_csv(path + 'GeneSymbol.csv',header=None,squeeze=True)
    df = (
        pd.read_csv(path + "parcelExpression.csv", header=None)
        .rename(columns={0:'label_id'})
        .set_index('label_id')
        .set_axis(genes, axis=1)
    )
    return df 

def read_labels_aurina(path, annotation):
    """
    Get labels from Aurina
    """
    labels = (
        pd.read_csv(path + 'SampleCoordinates.csv')
        .set_axis(['label', 'mni_x', 'mni_y', 'mni_z'], axis=1)
        .set_index(['mni_x', 'mni_y', 'mni_z'])
        .join(annotation.set_index(['mni_x', 'mni_y', 'mni_z']).loc[:,'well_id'])
        .set_index('well_id').sort_index()['label']
    )
    return labels
