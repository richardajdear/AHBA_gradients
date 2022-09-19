# Neurosynth metaanalysis with nimare

import numpy as np, pandas as pd
import nibabel as nib
from nilearn.maskers import NiftiLabelsMasker

from nimare.dataset import Dataset
from nimare.extract import fetch_neurosynth
from nimare.io import convert_neurosynth_to_dataset
from nimare import meta

import requests
import json
import re



def get_neurosynth_dset(dset_file = "../outputs/neurosynth_dset.pkl",
                        data_dir = '../data/neurosynth-data'):
    """
    Loads neurosynth dataset in nimare format from saved file.
    If saved dataset isn't found, attempts to download the neurosynth database 
    and convert into nimare format.
    """
    try:
        neurosynth_dset = Dataset.load(dset_file)
    except BaseException as err:
        print(f"Unexpected {err=}, {type(err)=}")
        print(f"Downloading neurosynth database to {data_dir}...")
        neurosynth_db = fetch_neurosynth(
            data_dir=data_dir, 
            source='abstract', 
            vocab='terms'
        )[0]

        print("Download succeeded. Converting...")
        neurosynth_dset = convert_neurosynth_to_dataset(
            coordinates_file=neurosynth_db["coordinates"],
            metadata_file=neurosynth_db["metadata"],
            annotations_files=neurosynth_db["features"],
        )
        
    return neurosynth_dset


def get_cogatlas_concepts(url=None):
    """ 
    Fetches list of concepts from the Cognitive Atlas
    """

    if url is None:
        url = 'https://cognitiveatlas.org/api/v-alpha/concept'

    req = requests.get(url)
    req.raise_for_status()
    concepts = set([f.get('name') for f in json.loads(req.content)])

    return concepts


def get_matching_terms(neurosynth_dataset, terms_list):
    """
    Gets list of terms with a match in neurosynth dataset labels
    """
    neurosynth_labels = neurosynth_dataset.get_labels()
    neurosynth_terms = [label.replace('terms_abstract_tfidf__', '') 
                        for label in neurosynth_labels]
    matched_terms = list(terms_list.intersection(neurosynth_terms))
    print(f"Matched {len(matched_terms)}/{len(terms_list)} \
            terms in Neurosynth vocab (total {len(neurosynth_terms)})")
    return matched_terms


def run_meta_one_term(term, neurosynth_dataset, estimator='MKDADensity'):
    """
    Runs meta analysis for one term
    """
    if estimator=='MKDADensity':
        meta_estimator = meta.cbma.mkda.MKDADensity()
    else:
        raise ValueError('Undefined meta estimator')

    label = f'terms_abstract_tfidf__{term}'
    label_positive_ids = neurosynth_dataset.get_studies_by_label(label, 0.05)
    label_negative_ids = list(set(neurosynth_dataset.ids) - set(label_positive_ids))

    if len(label_positive_ids)==0:
        raise ValueError('Term not matched')

    label_positive_dset = neurosynth_dataset.slice(label_positive_ids)
    label_negative_dset = neurosynth_dataset.slice(label_negative_ids)

    meta_result = meta_estimator.fit(label_positive_dset, label_negative_dset)

    return meta_result


def parcellate_meta(meta_result, 
                    parcellation_img_file="../data/parcellations/HCP-MMP_1mm.nii.gz"):
    """
    Parcellates a single meta analysis result with a specified parcellation image
    """

    zmap = meta_result.get_map('z')
    parcellation_img = nib.load(parcellation_img_file)
    parcellated_zmap = NiftiLabelsMasker(parcellation_img).fit_transform(zmap).squeeze()
    
    return parcellated_zmap


def get_parcellated_terms(terms, neurosynth_dataset, estimator='MKDADensity'):
    parcellated_maps = {}
    for i,term in enumerate(terms):
        try:
            print(f"({i}/{len(terms)}) Running meta analysis for '{term}'")
            _meta = run_meta_one_term(term, neurosynth_dataset, estimator)
        except ValueError:
            print('Skipping...')
            pass
        _meta_parcellated = parcellate_meta(_meta)
        parcellated_maps[term] = pd.Series(_meta_parcellated)

    return pd.concat(parcellated_maps, axis=1)


def list_disorder_terms():
    return [
        'anxiety',
        'disorder',
        'mdd','ptsd', 'ahdh', 'scz', 'sczd',
        'dementia',
        'depressive',
        'addiction',
        'psychopathology',
        'posttraumatic',
        'mood',
        'psychiatric',
        'bipolar',
        'compulsive',
        'disease',        
        'schizophreni',
        'psychosis',
        'parkinson',
        'stress',
        'degenerative',
        'autism'
    ]

def list_anatomical_terms():
    return [
        'cortex',
        'cortical',
        'cortico',
        'hypothalamus',
        'subgenual',
        'orbitofrontal',
        'amygdala',
        'paralimbic',
        'v1',
        'm1',
        'parieto',
        'parietal',
        'occipital',
        'pallidus',
        'middle',
        'globus',
        'regional',
        'gyrus',
        'sts',
        'corpus',
        'lobule',
        'ofc',
        'mtg',
        'callosum',
        'temporo'
        'mm',
        'early visual',
        'hippocampus',
        'accumbens',
        'primary visual',
        'striatum',
        'temporal',
        'ventral',
        'periaqueductal',
        'limbic',
        'medial',
        'anterior',
        'cortices',
        'v5',
        'ventromedial',
        'dorsolateral',
        'prefrontal',
        'pfc',
        'sma',
        'premotor',
        'supplementary',
        'dorsal',
        'caudal',
        'frontopolar',
        'lateral',
        'dlpfc',
        'dorsolateral',
        'primary motor',
        'superior',
        'precentral',
        'sensorimotor',
        'ipsilateral',
        'posterior',
        'contralateral',
        's1',
        'mtg',
        'somatosensory',
        'thalamus',
        'midbrain',
        'cingulate',
        'sts'
    ]