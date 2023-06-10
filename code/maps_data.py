# Functions to read and clean spatial data
import numpy as np, pandas as pd
from processing_helpers import *

def get_maps(data_dir="../data/cortical_maps.csv"):
    names = {
        'paquola_FC': 'Functional Hubs (FC)',

        'T1T2': 'Myelination (T1w/T2w)',
        'glasser_GI': 'Aerobic Glycolysis (AG)',
        'MEG_theta': 'MEG Theta (MEGθ)',
        'd_intramod60': 'Age 14-26 change in Short-Range\nIntramodular Connectivity (ΔShIm)',
        'dMT': 'Age 14-26 change in \n Myelination (ΔMT)',
    }

    maps = (pd.read_csv(data_dir)
            .drop(['d_intramod60','hill_dev','dCT'], axis=1)
            .apply(lambda x: (x-np.mean(x))/np.std(x))
            .set_index(get_labels_hcp()[:180])
    )

    return maps


def get_meg_maps(data_dir="../data/meg_HCPS1200.csv"):
    meg_maps = (pd.read_csv(data_dir, index_col=0)
                .drop(['myelinmap','thickness'], axis=1)
                .set_axis(get_labels_hcp()[:180])
              )
    return meg_maps


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


### LEGACY

# def get_ctmt_from_dk(project_to_hcp=True):
#     ctmt = (
#         pd.read_csv("../data/whitakervertes2016_complete.csv", index_col=0)
#         .query('hemi=="l"')
#         .set_index('label')
#         .loc[:, ['CT','CT_delta','MT','MT_delta']]#,'PLS2']]
#     )
    
#     if project_to_hcp:
#         ctmt = ctmt.pipe(dk_to_hcp).set_axis(get_labels_hcp()[:180])
            
#     return ctmt


# def get_brainchart_maps(data="../data/Peaks_Table_2_2.csv"):
#     maps = (
#         pd.read_csv(data, index_col=0)
#         .drop('Nboot',axis=1)
#         .set_index('feat')
#         .rename_axis('label')
#     )
#     maps = maps.set_index('lh_' + maps.index)
 
#     return maps


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