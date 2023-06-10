# Code for triplet analyses

from abagen import fetch_microarray
from itertools import combinations
from processing_helpers import *
from gradientVersion import *


def define_triplets(data_dir =  "../data/abagen-data/microarray/"):
    # Get donor ids
    files = fetch_microarray(donors='all', data_dir=data_dir)
    donors = list(files.keys())
    
    # Define triplets
    triplets_list = [list(x) for x in list(combinations(range(6), 3))]
    triplets_names = [''.join(map(str,x)) for x in triplets_list]
    triplets_dict = dict(zip(triplets_names, triplets_list))
    triplets_dict_donors = {k:[donors[i] for i in v] for k,v in triplets_dict.items()}
    disjoint_triplets = [list(x) for x in combinations(triplets_names,2) \
                         if not any(set(list(x[0])).intersection(set(list(x[1]))))]
    
    return triplets_dict_donors, disjoint_triplets


def get_triplets(**kwargs):
    # Get triplet definitions
    triplets_dict_donors, disjoint_triplets = define_triplets()
    
    triplets = {}
    for name, donors in triplets_dict_donors.items():
        expression, gene_stability = get_expression_abagen(
            donors = donors, 
            return_stability = True,
            verbose = 0,
            **kwargs
            )
        triplets[name] = gradientVersion().fit(expression)
        triplets[name].gene_stability = gene_stability
        print(f'Done triplet {name}')
    
    return triplets


def filter_triplet_stability(triplets, stability_threshold=0, **kwargs):
    # Iterate through each triplet
    triplets_filtered = {}

    for name, triplet in triplets.items():
        # Filter gene stability for the triplet
        gene_mask = triplet.gene_stability.rank(pct=True) > stability_threshold
        triplet_expression_filtered = triplet.expression.loc[:, gene_mask]
        
        # Don't use marker genes in case they are not retained in triplet
        triplets_filtered[name] = (gradientVersion(marker_genes=[], **kwargs)
                                .fit(triplet_expression_filtered, message=False))
        triplets_filtered[name].gene_stability = triplet.gene_stability
            
    return triplets_filtered


def get_triplets_stability_levels(triplets, **kwargs):
    """
    Compute gradients over different gene stability levels for each triplet
    """
    triplets_stability_levels = {}
    for stability in [i/10 for i in range(0,10)]:
         triplets_stability_levels[stability] = filter_triplet_stability(triplets, stability, **kwargs)
    print("Computed triplet stability")

    return triplets_stability_levels


def disjoint_corrs(triplets_version, how='coefs', match=True):
    """
    Correlate all pairs of disjoint triplets for a given parameter setting
    """
    # Get triplet definitions
    triplets_dict_donors, disjoint_triplets = define_triplets()
    
    corrs = {}
    for pair in disjoint_triplets:
#     for pair in all_pairs:
        name = '-'.join(pair)
        a = triplets_version[pair[0]]
        b = triplets_version[pair[1]]

        if how=='coefs':
            df_corr = a.corr_coefs(b, match=match)
        else:
            df_corr = a.corr_scores(b, match=match)
    
        if match:
            corrs[name] = df_corr['corr'].values
        else:
            corrs[name] = df_corr.pipe(np.diag)
            
    return pd.DataFrame(corrs)



def make_triplet_ds_plot(triplets_ds_levels, with_weights=False):
    """
    Compile the disjoint correlations from a set of triplets at different ds levels
    """
    
    scores_weights_dict = {}
    scores_dict = {name: disjoint_corrs(t, how='scores') for name, t in triplets_ds_levels.items()}
    scores_weights_dict['Region scores'] = pd.concat(scores_dict)
    
    if with_weights:
        weights_dict = {name: disjoint_corrs(t) for name, t in triplets_ds_levels.items()}
        scores_weights_dict['Gene weights'] = pd.concat(weights_dict)
        
    triplet_versions = (
        pd.concat(scores_weights_dict)
        .reset_index()
        .rename(columns={'level_0':'how', 'level_1':'version', 'level_2':'component'})
        .assign(component=lambda x: x['component'].replace({i:f'{i+1}' for i in range(5)}))
        .melt(id_vars=['how', 'version', 'component'], var_name='pair', value_name='corr')
        .assign(corr_abs = lambda x: np.abs(x['corr']))
    )
    
    return triplet_versions


def make_triplet_versions_plot(triplet_versions_dict, with_weights=False):    
    """
    xx
    """
    scores_weights_dict = {}
    scores_dict = {name: disjoint_corrs(version, how='scores') for name, version in triplet_versions_dict.items()}
    scores_weights_dict['Region scores'] = pd.concat(scores_dict)
    
    # if with_weights:
    #     weights_dict = {name: disjoint_corrs(version) for name, version in triplet_versions_dict.items()}
    #     scores_weights_dict['Gene weights'] = pd.concat(weights_dict)

    triplet_versions = (
        pd.concat(scores_weights_dict)
        .reset_index()
        .rename(columns={'level_0':'how', 'level_1':'version', 'level_2':'component'})
        .assign(component=lambda x: x['component'].replace({i:f'{i+1}' for i in range(5)}))
        .melt(id_vars=['how', 'version', 'component'], var_name='pair', value_name='corr')
        .assign(corr_abs = lambda x: np.abs(x['corr']))
        .assign(version = lambda x: pd.Categorical(x['version'], categories=x['version'].unique(),ordered=True))
    )
    
    return triplet_versions



# def filter_triplet_ds(triplets, ds_threshold=0, **kwargs):

#     # Iterate through each triplet
#     triplets_ds = {}
#     for name, triplet in triplets.items():
        
#         # Filter gene stability for the triplet
#         mask = triplet.stability.rank(pct=True) > ds_threshold
#         triplet_expression_ds = triplet.expression.loc[:, mask]
        
#         # Don't use marker genes in case they are not retained in triplet
#         triplets_ds[name] = (gradientVersion(marker_genes=[], **kwargs)
#                                 .fit(triplet_expression_ds, message=False))
            
#     return triplets_ds

# def get_triplets_ds_levels(triplets, **kwargs):
#     """
#     x
#     """
#     triplets_ds_levels = {}
#     for ds in [i/10 for i in range(0,10)]:
#          triplets_ds_levels[ds] = filter_triplet_ds(triplets, ds, **kwargs)

#     return triplets_ds_levels