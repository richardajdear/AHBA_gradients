# Functions to read and clean spatial data
import numpy as np, pandas as pd
from processing import *

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

def get_yeo_mesulam(scores=None):
    lobe_colors = {'Occ':'#E41A1C', 'Fr':'#4A72A6', 'Par':'#48A462', 'Temp':'#7E6E85', 'Ins':'#D16948'}
    lobe_names = {'Occ':'Occipital', 'Fr':'Frontal', 'Par':'Parietal', 'Temp':'Temporal', 'Ins':'Insula'}
    
    mesulam_colors = {1:"#FFFFB3", 2:"#FB8072", 3:"#BEBADA", 4:"#80B1D3"}
    mesulam_names = {1:'Paralimbic',2:'Heteromodal',3:'Unimodal',4:'Idiotypic'}

    hcp_yeo_mesulam = (pd.read_csv('../data/parcellations/HCP-MMP1_UniqueRegionList.txt')
     .assign(id=lambda x: [id-20 if id>180 else id for id in x['regionID']])
     .set_index('id')
     .loc[:180, ['region', 'Lobe']]
     .assign(
        Lobe_colors = lambda x: x['Lobe'].map(lobe_colors),
        Lobe_names = lambda x: x['Lobe'].map(lobe_names)
     )
     .join(pd.read_csv("../data/yeo_asymmetric.csv").set_index('id').drop('label', axis=1))
     .join(pd.read_csv("../data/mesulam_hcp.csv").set_index('id').drop('label', axis=1))
     .assign(
        Mesulam_colors = lambda x: x['Mesulam'].map(mesulam_colors),
        Mesulam_names = lambda x: x['Mesulam'].map(mesulam_names)
     )
    )


    if scores is not None:
        hcp_yeo_mesulam = hcp_yeo_mesulam.join(scores)
    
    return hcp_yeo_mesulam

    # hcp_yeo_mesulam.to_csv("../data/hcp_yeo_mesulam.csv")


# def get_mesulam_ve_yeo(data_dir="../data/mesulam_ve_yeo.csv"):
#     short_names_dict = {
#         'Mesulam_names': {
#             'Heteromodal':'Het.',
#             'Unimodal':'Uni.',
#             'Paralimbic':'Par.',
#             'Idiotypic':'Idi.'
#         },
#         'Yeo_names': {
#             'Visual':'VIS',
#             'Somatomotor':'SMN',
#             'Dorsal attention':'DAN',
#             'Central attention':'CAN',
#             'Default mode':'DMN',
#             'Frontoparietal':'FPN',
#             'Limbic':'LIM'
#         },
#     }

#     mesulam_ve_yeo = (pd.read_csv(data_dir, index_col=0)
#                       .drop(['x','y','z'], axis=1)
#                       .replace({'label':{'7PL':'7Pl'}})
#                       .replace({'Mesulam':{4:1,2:4,3:2}})
#                     #   .replace({'Yeo_names':{'Central attention':'Central atten.','Dorsal attention':'Dorsal atten.'}})
#                       .replace(short_names_dict)
#                       .set_index('label'))
#     return mesulam_ve_yeo


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