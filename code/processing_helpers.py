# Helper functions for processing

import numpy as np, pandas as pd
import abagen
from abagen.correct import keep_stable_genes
import abagen_allen_tweaked
import nibabel as nib
from pcaVersion import pcaVersion


def get_expression_abagen(atlas, save_name=None, verbose=0, 
                          only_left=True,
                          only_cortex=True,
                          ibf_threshold=0.5,
                          probe_selection='diff_stability', 
                          lr_mirror='rightleft',
                          region_agg='donors',
                          tolerance=2,
                          sample_norm='srs',
                          gene_norm='srs',
                          DS_threshold=0,
                          return_donors=False,
                          return_counts=False,
                          return_labels=False,
                          return_stability=False,
                          **kwargs):
    """
    Get expression matrix and matching labels from Abagen X
    """

    out = abagen_allen_tweaked.get_expression_data(
        atlas=atlas['image'],
        atlas_info=atlas['info'],
        verbose=verbose,
        only_left=only_left,
        only_cortex=only_cortex,
        ibf_threshold=ibf_threshold,
        probe_selection=probe_selection,
        lr_mirror=lr_mirror,
        region_agg=region_agg,
        tolerance=tolerance,
        sample_norm=sample_norm,
        gene_norm=gene_norm,
        data_dir='../data/abagen-data/microarray',
        n_proc=32,
        return_counts=return_counts,
        return_labels=return_labels,
        return_donors=True, # always true to get DS
        **kwargs
    )
    
    if return_labels or return_counts:
        expression_all_genes = out[0]
    else:
        expression_all_genes = out
    
    expression, stability = keep_stable_genes(expression_all_genes, threshold=DS_threshold, return_stability=True)
    stability = pd.Series(stability, index=expression_all_genes[0].columns)
    print(f'{expression[0].shape[1]} genes remain after filtering for top {round(1-DS_threshold,2)} differential stability')
    
    if not return_donors:
        expression = pd.concat(expression).groupby(level=0).mean()
    
    if save_name is not None:
        expression.to_csv("../data/abagen-data/expression/" + save_name + ".csv")
        if return_labels:
            out[1].to_csv("../data/abagen-data/labels/" + save_name + ".csv")
    
    out_ = (expression,)
    if return_labels:
        out_ += (out[1],)
    if return_counts:
        if return_labels:
            out_ += (out[2],)
        else:
            out_ += (out[1],)
    if return_stability:
        out_ += (stability,)
    if len(out_) == 1:
        out_ = out_[0]
    
    return out_

    
    

def fetch_dk(native=True, only_cortex=False, only_left=False):
    """
    Get DK atlas and remove subcortex
    """
    atlas_dk = abagen.fetch_desikan_killiany(native=native)
    
    def remove_subcortex(image):
        img = nib.load(image)
        if only_left:
            img_cortex = np.where((img.dataobj[:] > 34), 0, img.dataobj[:])
        else:
            img_cortex = np.where((img.dataobj[:] > 34) & (img.dataobj[:] < 48), 0, img.dataobj[:])
        return img.__class__(img_cortex, img.affine, img.header)
        
    if only_cortex:
        for donor, image in atlas_dk['image'].items():
            atlas_dk['image'][donor] = remove_subcortex(image)
#         atlas_dk['info'] = pd.read_csv(atlas_dk['info']).replace('subcortex/brainstem','other')
    
    return atlas_dk



def fetch_hcp(native=True, only_left=False):
    """
    Get HCP atlas
    """
    hcp_info = (
        pd.read_csv('../data/HCP-MMP1_UniqueRegionList.txt')
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
            'image':'../data/HCP-MMP_1mm.nii.gz',
            'info':pd.read_csv('../data/HCP-MMP1_UniqueRegionList.txt').rename(columns={'LR':'hemisphere', 'regionID':'id', 'region':'label'}).assign(structure='cortex')
        }

    def remove_right(image):
        img = nib.load(image)
        img_cortex = np.where((img.dataobj[:] > 180), 0, img.dataobj[:])
        return img.__class__(img_cortex, img.affine, img.header)
    
    if only_left:
        for donor, image in atlas_hcp['image'].items():
            atlas_hcp['image'][donor] = remove_right(image)

    return atlas_hcp



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
        pd.read_csv('../data/HCP-MMP1_UniqueRegionList.txt')
        .rename(columns={'LR':'hemisphere', 'regionID':'id', 'region':'label'})
        .assign(structure='cortex')
        .assign(id=lambda x: [id-20 if id>180 else id for id in x['id']]) 
        # Aurina's HCP images code regions as 1-360, not 1-180,201-380
        # so recode to match the image files
    )
    labels_hcp = hcp_info.set_index('id')['label']
    return labels_hcp

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
