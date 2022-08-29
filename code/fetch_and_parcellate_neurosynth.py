# -*- coding: utf-8 -*-
"""
Copied from https://github.com/netneurolab/ipn-summer-school/blob/main/lectures/2021-07-02/13-15/fetch_and_parcellate_neurosynth.py
and edited


Script for performing NeuroSynth-style meta-analyses for all available
Cognitive Atlas concepts. Script modified from https://github.com/netneurolab/
markello_spatialnulls/blob/master/scripts/empirical/fetch_neurosynth_maps.py,
for the IPN Summer School > Quantitative and Computational Neuroscience >
Advanced Analytics for Neuroscience > Contextualizing Results lecture.

Again, this is mostly Ross' code.
"""

import contextlib
import json
import os
from pathlib import Path
import warnings

import pandas as pd
import numpy as np
import requests
from nilearn.input_data import NiftiLabelsMasker
from nilearn._utils import check_niimg

# neurosynth downloaded as
# pip install git+https://github.com/neurosynth/neurosynth
import neurosynth as ns
from nilearn.datasets import fetch_atlas_schaefer_2018


# /sigh
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# this is where the raw and parcellated data will be stored
NSDIR = Path('../data/neurosynth-data/raw').resolve()
PARDIR = Path('../data/neurosynth-data/parcellated').resolve()

# these are the images from the neurosynth analyses we'll save
# can add 'uniformity-test_z' plus more, if desired
IMAGES = ['association-test_z']


def fetch_ns_data(directory):
    """ Fetches NeuroSynth database + features to `directory`
    Paramerters
    -----------
    directory : str or os.PathLike
        Path to directory where data should be saved
    Returns
    -------
    database, features : PathLike
        Paths to downloaded NS data
    """

    directory = Path(directory)

    # if not already downloaded, download the NS data and unpack it
    database, features = directory / 'database.txt', directory / 'features.txt'
    if not database.exists() or not features.exists():
        with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):
            ns.dataset.download(path=directory, unpack=True)
        try:  # remove tarball if it wasn't removed for some reason
            (directory / 'current_data.tar.gz').unlink()
        except FileNotFoundError:
            pass

    return database, features


def get_cogatlas_concepts(url=None):
    """ Fetches list of concepts from the Cognitive Atlas
    Parameters
    ----------
    url : str
        URL to Cognitive Atlas API
    Returns
    -------
    concepts : set
        Unordered set of terms
    """

    if url is None:
        url = 'https://cognitiveatlas.org/api/v-alpha/concept'

    req = requests.get(url)
    req.raise_for_status()
    concepts = set([f.get('name') for f in json.loads(req.content)])

    return concepts


def run_meta_analyses(database, features, use_features=None, outdir=None):
    """
    Runs NS-style meta-analysis based on `database` and `features`
    Parameters
    ----------
    database, features : str or os.PathLike
        Path to NS-style database.txt and features.txt files
    use_features : list, optional
        List of features on which to run NS meta-analyses; if not supplied all
        terms in `features` will be used
    outdir : str or os.PathLike
        Path to output directory where derived files should be saved
    Returns
    -------
    generated : list of str
        List of filepaths to generated term meta-analysis directories
    """

    # check outdir
    if outdir is None:
        outdir = NSDIR
    outdir = Path(outdir)

    # make database and load feature names; annoyingly slow
    dataset = ns.Dataset(str(database))
    dataset.add_features(str(features))
    features = set(dataset.get_feature_names())

    # if we only want a subset of the features take the set intersection
    if use_features is not None:
        features = set(features) & set(use_features)
    pad = max([len(f) for f in features])

    generated = []
    for word in sorted(features):
        msg = f'Running meta-analysis for term: {word:<{pad}}'
        print(msg, end='\r', flush=True)

        # run meta-analysis + save specified outputs (only if they don't exist)
        path = outdir / word.replace(' ', '_')
        path.mkdir(exist_ok=True)
        if not all((path / f'{f}.nii.gz').exists() for f in IMAGES):
            ma = ns.MetaAnalysis(dataset, dataset.get_studies(features=word))
            ma.save_results(path, image_list=IMAGES)

        # store MA path
        generated.append(path)

    print(' ' * len(msg) + '\b' * len(msg), end='', flush=True)

    return generated


def parcellate_meta(outputs, annots, fname, regions):
    # empty dataframe to hold our parcellated data
    data = pd.DataFrame(index=regions)
    mask = NiftiLabelsMasker(annots, resampling_target='data')

    for outdir in outputs:
        cdata = []
        mgh = outdir / 'association-test_z.nii.gz'

        cdata.append(mask.fit_transform(
            check_niimg(mgh.__str__(), atleast_4d=True)).squeeze())

        # aaaand store it in the dataframe
        data = data.assign(**{outdir.name: np.hstack(cdata)})

    # now we save the dataframe! wooo data!
    data.to_csv(fname, sep=',')
    return fname


if __name__ == '__main__':
    NSDIR.mkdir(parents=True, exist_ok=True)
    PARDIR.mkdir(parents=True, exist_ok=True)

    # get concepts from CogAtlas and run relevant NS meta-analysess,
    database, features = fetch_ns_data(NSDIR)
    generated = run_meta_analyses(database, features, get_cogatlas_concepts(),
                                  outdir=NSDIR)

    # get parcellations that we'll use to parcellate data
    schaefer = fetch_atlas_schaefer_2018(n_rois=200, resolution_mm=2)
    labels = []
    for i in range(len(schaefer['labels'])):
        labels.append(schaefer['labels'][i].decode("utf-8"))

    # parcellate data and save to directory
    parcellate_meta(generated, schaefer['maps'],
                    PARDIR / 'atl-schaefer2018_res-200_neurosynth.csv',
                    regions=labels)