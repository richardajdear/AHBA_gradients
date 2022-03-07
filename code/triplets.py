# Code for triplet analyses

from abagen import fetch_microarray
from itertools import combinations
from processing_helpers import *
from pcaVersion import *
from gradientVersion import *


def define_triplets(data_dir =  "~/rds/rds-cam-psych-transc-Pb9UGUlrwWc/Cam_LIBD/AHBA_data/abagen-data/microarray/"):
    # Get donor ids
    files = fetch_microarray(donors='all', data_dir=data_dir)
    donors = list(files.keys())
    
    # Define triplets
    triplets_list = [list(x) for x in list(combinations(range(6), 3))]
    triplets_names = [''.join(map(str,x)) for x in triplets_list]
    triplets_dict = dict(zip(triplets_names, triplets_list))
    triplets_dict_donors = {k:[donors[i] for i in v] for k,v in triplets_dict.items()}
    disjoint_triplets = [list(x) for x in combinations(triplets_names,2) if not any(set(list(x[0])).intersection(set(list(x[1]))))]
    
    return triplets_dict_donors, disjoint_triplets


def get_triplets(**kwargs):
    # Get triplet definitions
    triplets_dict_donors, disjoint_triplets = define_triplets()
    
    triplets = {}
    for name, donors in triplets_dict_donors.items():
        expression, stability = get_expression_abagen(donors=donors, **kwargs, return_stability=True)
        triplets[name] = pcaVersion(expression, message=False, scale=False) # don't scale so that we can use gradientVersion
        triplets[name].stability = stability
        print(f'Done triplet {name}')
    
    return triplets


def filter_triplet_ds(triplets, ds_threshold=0, only_brain=False, use_gradientVersion=False, **kwargs):
    # Get brain genes
    if only_brain:
        brain_genes = pd.read_csv("../data/jakob_brain_genes.csv", squeeze=True)
    
    # Iterate through each triplet
    triplets_ds = {}
    for name, triplet in triplets.items():
        
        # Filter DS for the triplet
        mask = triplet.stability.rank(pct=True) > ds_threshold
        triplet_expression_ds = triplet.expression.loc[:, mask]
        
        # Filter for only brain genes
        if only_brain:
            brain_mask = set(brain_genes).intersection(set(triplet_expression_ds.columns))
            triplet_expression_ds = triplet_expression_ds.loc[:, brain_mask]
        
        # Use pcaVersion, or use gradientMaps
        if use_gradientVersion:
            triplets_ds[name] = gradientVersion(**kwargs).fit(triplet_expression_ds, message=False)
        else:
            triplets_ds[name] = pcaVersion(triplet_expression_ds, message=False, **kwargs)
            
    return triplets_ds


def disjoint_corrs(triplets_version, how='coefs', match=True):
    # Get triplet definitions
    triplets_dict_donors, disjoint_triplets = define_triplets()
    
    corrs = {}
    for pair in disjoint_triplets:
#     for pair in all_pairs:
        name = '-'.join(pair)
        pca1 = triplets_version[pair[0]]
        pca2 = triplets_version[pair[1]]
    
        if how=='coefs':
            df_corr = pca1.corr_coefs(pca2, match=match)
        else:
            df_corr = pca1.corr_scores(pca2, match=match)
    
        if match:
            corrs[name] = df_corr['corr'].values
        else:
            corrs[name] = df_corr.pipe(np.diag)
            
    return pd.DataFrame(corrs)


def get_triplets_ds_levels(triplets, **kwargs):
    """
    x
    """
    triplets_ds_levels = {}
    for ds in [i/10 for i in range(0,10)]:
         triplets_ds_levels[ds] = filter_triplet_ds(triplets, ds, **kwargs)

    return triplets_ds_levels



def make_triplet_versions_plot(triplet_versions_dict, with_coefs=True):    
    
    scores_coefs_dict = {}
    scores_dict = {name: disjoint_corrs(version, how='scores') for name, version in triplet_versions_dict.items()}
    scores_coefs_dict['Region scores'] = pd.concat(scores_dict)
    
    if with_coefs:
        coefs_dict = {name: disjoint_corrs(version) for name, version in triplet_versions_dict.items()}
        scores_coefs_dict['Gene weights'] = pd.concat(coefs_dict)

    triplet_versions = (
        pd.concat(scores_coefs_dict)
        .reset_index()
        .rename(columns={'level_0':'how', 'level_1':'version', 'level_2':'component'})
        .assign(component=lambda x: x['component'].replace({i:f'{i+1}' for i in range(5)}))
        .melt(id_vars=['how', 'version', 'component'], var_name='pair', value_name='corr')
        .assign(corr_abs = lambda x: np.abs(x['corr']))
        .assign(version = lambda x: pd.Categorical(x['version'], categories=x['version'].unique(),ordered=True))
    )
    
    return triplet_versions


def make_triplet_ds_plot(triplets_ds_levels, with_coefs=True):
    """
    x
    """
    
    scores_coefs_dict = {}
    scores_dict = {name: disjoint_corrs(t, how='scores') for name, t in triplets_ds_levels.items()}
    scores_coefs_dict['Region scores'] = pd.concat(scores_dict)
    
    if with_coefs:
        coefs_dict = {name: disjoint_corrs(t) for name, t in triplets_ds_levels.items()}
        scores_coefs_dict['Gene weights'] = pd.concat(coefs_dict)
        
    triplet_versions = (
        pd.concat(scores_coefs_dict)
        .reset_index()
        .rename(columns={'level_0':'how', 'level_1':'version', 'level_2':'component'})
        .assign(component=lambda x: x['component'].replace({i:f'{i+1}' for i in range(5)}))
        .melt(id_vars=['how', 'version', 'component'], var_name='pair', value_name='corr')
        .assign(corr_abs = lambda x: np.abs(x['corr']))
    )
    
    return triplet_versions