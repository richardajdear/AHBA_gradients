# Helper functions for processing brainspan data

import numpy as np, pandas as pd
from processing_helpers import *

def get_brainspan(bs_dir = "../data/brainspan-data/gene_matrix_rnaseq/"):
    """
    Read in brainspan data as expression matrix, column labels, and row labels
    """
    # Expression matrix
    bs_exp = pd.read_csv(bs_dir + "expression_matrix.csv", header=None, index_col=0)
    # Columns
    bs_col = (pd.read_csv(bs_dir + "columns_metadata.csv")
              .assign(age = lambda x: pd.Categorical(x['age'], ordered=True, categories=x['age'].unique()))
             )
    # Rows
    bs_row = pd.read_csv(bs_dir + "rows_metadata.csv")
    return bs_exp, bs_col, bs_row


def get_age_groups():
    """
    Define age groupings for use in Brainspan analysis
    """
    age_groups = {        
        # 'Best' groupings ...
        # '12 pcw': '12-17 pcw',
        # '13 pcw': '12-17 pcw',
        # '16 pcw': '12-17 pcw',
        # '17 pcw': '12-17 pcw',
        # '19 pcw': '19-37 pcw',
        # '21 pcw': '19-37 pcw',
        # '24 pcw': '19-37 pcw',
        # '37 pcw': '19-37 pcw',
        # '4 mos': '4 mos',
        # '1 yrs': '1 yrs',
        # '2 yrs': '2 yrs',
        # '3 yrs': '3 yrs',
        # '8 yrs': '8-13 yrs',
        # '11 yrs': '8-13 yrs',
        # '13 yrs': '8-13 yrs',
        # '18 yrs': '18-40 yrs',
        # '19 yrs': '18-40 yrs',
        # '21 yrs': '18-40 yrs',
        # '23 yrs': '18-40 yrs',
        # '30 yrs': '18-40 yrs',
        # '36 yrs': '18-40 yrs',
        # '37 yrs': '18-40 yrs',
        # '40 yrs': '18-40 yrs',   
        
        # # Simple groupings ...
        '8 pcw': 'Pre-Birth',
        '12 pcw': 'Pre-Birth',
        '13 pcw': 'Pre-Birth',
        '16 pcw': 'Pre-Birth',
        '17 pcw': 'Pre-Birth',
        '19 pcw': 'Pre-Birth',
        '21 pcw': 'Pre-Birth',
        '24 pcw': 'Pre-Birth',
        '37 pcw': 'Pre-Birth',
        '4 mos': 'Birth-13 yrs',
        '10 mos': 'Birth-13 yrs',
        '1 yrs': 'Birth-13 yrs',
        '2 yrs': 'Birth-13 yrs',
        '3 yrs': 'Birth-13 yrs',
        '4 yrs': 'Birth-13 yrs',
        '8 yrs': 'Birth-13 yrs',
        '11 yrs': 'Birth-13 yrs',
        '13 yrs': 'Birth-13 yrs',
        '18 yrs': '18-40 yrs',
        '19 yrs': '18-40 yrs',
        '21 yrs': '18-40 yrs',
        '23 yrs': '18-40 yrs',
        '30 yrs': '18-40 yrs',
        '36 yrs': '18-40 yrs',
        '37 yrs': '18-40 yrs',
        '40 yrs': '18-40 yrs',   
    }
    return age_groups


def get_bs_cortex_mapping():
    """
    Define mapping between Brainspan regions and HCP cortex groups
    """
    bs_cortex_mapping = {
        'primary visual cortex (striate cortex, area V1/17)': 'Primary_Visual',
        'posteroventral (inferior) parietal cortex': 'Inferior_Parietal',
        'primary somatosensory cortex (area S1, areas 3,1,2)': 'Somatosensory',
        'primary motor cortex (area M1, area 4)': 'Motor',
        'dorsolateral prefrontal cortex': 'Dorsolateral_Prefrontal',
        'ventrolateral prefrontal cortex': 'Inferior_Frontal',    
        'anterior (rostral) cingulate (medial prefrontal) cortex': 'Anterior_Cingulate_and_Medial_Prefrontal',
        'orbital frontal cortex': 'Orbital_and_Polar_Frontal',
        'inferolateral temporal cortex (area TEv, area 20)': 'Lateral_Temporal',
        'primary auditory cortex (core)': 'Early_Auditory',
        'posterior (caudal) superior temporal cortex (area 22c)': 'Auditory_Association'
    }
    
    return bs_cortex_mapping


def get_hcp_bs_mapping_v2():
    """
    Define mapping between individual HCP regions and BrainSpan
    """
    hcp_bs_mapping = (
        pd.read_csv("../data/hcp_bs_mapping_v2.csv", index_col=None)
        # .query("keep==1")
        .assign(structure_name = lambda x: np.where(x['keep']==1, x['structure_name'], np.nan))
               )
    
    return hcp_bs_mapping

def get_hcp_bs_mapping_v3():
    """
    Define mapping between individual HCP regions and BrainSpan
    """
    hcp_bs_mapping = (
        pd.read_csv("../data/hcp_bs_mapping_v3.csv", index_col=None)
        # .query("keep==1")
        .assign(structure_name = lambda x: np.where(x['keep']==1, x['structure_name'], np.nan))
               )
    
    return hcp_bs_mapping


def get_dk_bs_mapping():
    """
    Define mapping between Brainspan regions and HCP cortex groups
    """
    dk_bs_mapping = {
        'pericalcarine':'primary visual cortex (striate cortex, area V1/17)',
        'inferior parietal':'posteroventral (inferior) parietal cortex',
        'postcentral':'primary somatosensory cortex (area S1, areas 3,1,2)',
        'precentral':'primary motor cortex (area M1, area 4)',
        'caudal middle frontal':'dorsolateral prefrontal cortex',
        'rostral middle frontal':'dorsolateral prefrontal cortex',
        'pars orbitalis':'ventrolateral prefrontal cortex',
        'pars opercularis':'ventrolateral prefrontal cortex',
        'pars triangularis':'ventrolateral prefrontal cortex',
        'rostral anterior cingulate':'anterior (rostral) cingulate (medial prefrontal) cortex',
        'caudal anterior cingulate':'anterior (rostral) cingulate (medial prefrontal) cortex',
        'frontal pole':'orbital frontal cortex',
        'lateral orbitofrontal':'orbital frontal cortex',
        'medial orbitofrontal':'orbital frontal cortex',
        'inferior temporal':'inferolateral temporal cortex (area TEv, area 20)',
        'fusiform':'inferolateral temporal cortex (area TEv, area 20)',
        'transverse temporal':'primary auditory cortex (core)',
        'superior temporal':'posterior (caudal) superior temporal cortex (area 22c)'
    }

    dk_bs_mapping = (
        pd.DataFrame.from_dict(dk_bs_mapping, orient='index')
        .assign(label = lambda x: 'lh_' + x.index.str.replace(' ',''))
        .set_index('label')
        .reindex(get_labels_dk()[:34])
        .set_axis(['structure_name'], axis=1)
        .assign(structure_name = lambda x: x['structure_name'].astype('string'))
        .reset_index()
    )
    
    return dk_bs_mapping



def get_short_bs_names():
    """
    Short Brainspan region names for plotting
    """
    bs_structure_name_short = {
        'primary visual cortex (striate cortex, area V1/17)': 'primary visual cortex',
        'posteroventral (inferior) parietal cortex': 'inferior parietal cortex',
        'primary somatosensory cortex (area S1, areas 3,1,2)': 'primary somatosensory cortex',
        'primary motor cortex (area M1, area 4)': 'primary motor cortex',
        'dorsolateral prefrontal cortex': 'dorsolateral prefrontal cortex',
        'ventrolateral prefrontal cortex': 'ventrolateral prefrontal cortex',
        'anterior (rostral) cingulate (medial prefrontal) cortex': 'medial prefrontal cortex',
        'orbital frontal cortex': 'orbital frontal cortex',
        'inferolateral temporal cortex (area TEv, area 20)': 'inferolateral temporal cortex',
        'primary auditory cortex (core)': 'primary auditory cortex',
        'posterior (caudal) superior temporal cortex (area 22c)':'posterior superior temporal cortex'
    }
    return bs_structure_name_short


def get_hcp_bs_mapping(hcp_info_file = "../data/parcellations/HCP-MMP1_UniqueRegionList.txt"):
    """
    Get HCP regions mapped to Brainspan from Brainspan cortex mapping
    """
    # Read HCP info file
    hcp_info = pd.read_csv(hcp_info_file)
    # Split somatosensory-Motor cortex
    hcp_info = (hcp_info.assign(cortex = lambda x: np.select(
        condlist = [(x['cortex'] == 'Somatosensory_and_Motor') & (x['region'] == '4'), 
                    (x['cortex'] == 'Somatosensory_and_Motor') & (x['region'] != '4')],
        choicelist = ['Motor', 'Somatosensory'], 
        default = x['cortex']
    )))
    # Get bs cortex mapping and invert
    bs_cortex_mapping = get_bs_cortex_mapping()
    cortex_bs_mapping = {v:k for k,v in bs_cortex_mapping.items()}
    # Explode BS regions to all HCP regions
    hcp_bs_mapping = (hcp_info
     .loc[lambda x: x['LR'] == 'L', ['region', 'cortex']]
     .assign(structure_name = lambda x: x['cortex'].map(cortex_bs_mapping).astype('string'))
     .assign(structure_name_short = lambda x: x['structure_name'].map(get_short_bs_names()))
     .sort_values('structure_name')
     .rename({'region':'label'}, axis=1)
                     )
    return hcp_bs_mapping
    

def get_mapped_scores(version, version_to_bs_mapping, mean=True):
    """
    Get gradient scores in mapped Brainspan regions
    """
    # Get gradients filtered to HCP regions matched in brainspan
    scores_filtered = (
     version.clean_scores()
     .set_index('label')
     .join(version_to_bs_mapping.set_index('label')['structure_name'])
     .dropna(axis=0)
    )

    scores_mean = (scores_filtered
     .groupby('structure_name')
     .mean()
     .apply(lambda x: (x-np.mean(x))/np.std(x))
    )
    
    if mean:
        return scores_mean
    else:
        return scores_filtered

    
# def get_mapped_pcs(pc_version, hcp_base, hcp_bs_mapping):
#     """
#     Get PC scores in mapped Brainspan regions
#     """
#     # Get PCs filtered to HCP regions matched in brainspan
#     pcs_filtered = (
#      hcp_base.score_from(pc_version)
#      .join(get_labels_hcp())
#      .rename_axis('id')
#      .join(hcp_bs_mapping.set_index('region'), on='label')
#      .dropna(axis=0)
#     )

#     pcs_cortex = (pcs_filtered
#      .groupby('cortex')
#      .mean()
#      .apply(lambda x: (x-np.mean(x))/np.std(x))
#     )
    
#     return pcs_filtered, pcs_cortex


def count_samples(bs_col, bs_cortex_mapping):
    """
    Count samples from each Brainspan region and age
    """
    age_region_counts = (bs_col
                         .assign(cortex = lambda x: x['structure_name'].map(bs_cortex_mapping))
     .pivot_table(values='column_num', index='cortex', columns='age', aggfunc='count')
                        )
    return age_region_counts



def clean_brainspan(bs_exp, bs_col, bs_row, bs_mapping):
    """
    Clean up Brainspan data into dataframe
    """
    mapped_regions = bs_mapping['structure_name'].dropna().unique()
    
    # Join region data
    # Filter for only mapped regions (i.e. no subcortex)
    bs_mapped = (pd.concat([
        bs_exp.T,
        bs_col.set_index('column_num')[['donor_id', 'age', 'structure_name']]
    ], axis=1)
     .loc[lambda x: np.isin(x['structure_name'], mapped_regions)] # filter for mapped regions
     .set_index(['donor_id', 'age', 'structure_name'])
     # .dropna(how='all')
    )
    
    # Join gene data and clean NAs and duplicates
    bs_clean = (bs_mapped
     .set_axis(bs_row['gene_symbol'], axis=1)
     .fillna(0)
     .loc[:, lambda x: (x != 0).all(axis=0)] # drop zero and na columns
     .loc[:, lambda x: ~x.columns.duplicated()] # some columns are duplicates
    )
    
    return bs_clean


def aggregate_brainspan_by_age(bs_clean, normalize=True):
    """
    Aggregate Brainspan donor brains by ages
    Drop age groups with <=3 samples
    """
    # Drop brains with <= 3 missing samples
    bs_donor_counts = bs_clean.groupby(['donor_id']).size()
    bs_keep = bs_donor_counts > 3
    ix_keep = pd.IndexSlice[bs_keep[bs_keep].index.get_level_values('donor_id'), :]
    bs_clean = bs_clean.loc[ix_keep,:]

    # Optionally normalize by donor
    if normalize:
        bs_clean = (bs_clean
         .groupby(['donor_id', 'age'])
         .apply(lambda x: (x-np.mean(x))/np.std(x))
        )
        
    bs_agg = (bs_clean
      .groupby(['age', 'structure_name'], observed=True).mean())
    return bs_agg



def compute_brainspan_scores(bs_agg, version, normalize=True):
    """
    Compute scores of AHBA gradients on donor-aggregated Brainspan
    Add in missing NA rows
    Optionally normalize by age
    """
    # Find matching genes
    gene_mask = version.weights.index.intersection(bs_agg.columns)
    
    # Score PCs
    bs_scores = (bs_agg.loc[:, gene_mask] @ version.weights.loc[gene_mask, :]).iloc[:, :3].set_axis(['G1','G2','G3'], axis=1)
    
    # Add missing regions to each age as NA
    bs_regions_index = pd.CategoricalIndex(
        bs_scores.index.get_level_values(1).unique().dropna()
    )
    bs_scores = (bs_scores
                 .groupby('age')
                 .apply(lambda x: x.droplevel(0).reindex(bs_regions_index))
                 )
    
    # Normalize by age
    if normalize:
        bs_scores = bs_scores.groupby('age').apply(lambda x: (x-np.mean(x))/np.std(x))
        
    return bs_scores


def correlate_brainspan_scores(bs_scores, ahba_scores, age_groups=None, rolling=None, plot=True):
    """
    Correlate Brainspan gradient scores with HCP gradient scores
    Take absolute correlation because PLS may invert scores
    Either by all ages, or in defined age groups
    Optionally melt for plotting
    """
    # Aggregate into age groups if desired
    if age_groups is not None:
        bs_scores = (bs_scores
        .reset_index()
        .assign(age_group = lambda x: x['age'].map(age_groups))
        .assign(age = lambda x: pd.Categorical(x['age_group'], ordered=True, categories = x['age_group'].unique()))
        .drop('age_group', axis=1)
        .groupby(['age', 'structure_name']).mean() # agg into age groups
                 )

    # Rolling average over ages if desired
    if rolling is not None:
        bs_scores = (bs_scores
        .reset_index()
        .groupby(['structure_name'])
        .rolling(rolling, center=False, on='age', min_periods=1).mean()
        .set_index(['age'], append=True)
        .droplevel(level=1)
                 )
        
    # Correlate
    # Take absolute because PLS may invert
    bs_scores_corr = (bs_scores
    .groupby('age').corrwith(ahba_scores)
    .abs() # take absolute
                  ) 
    # Clean up for plotting
    if plot:
        bs_scores_corr = (bs_scores_corr
                       .melt(var_name='G', value_name='corr', ignore_index=False)
                       .reset_index())
    return bs_scores_corr


def combine_scores(bs_scores, ahba_scores_mapped, age_groups):
    """
    Combine AHBA and Brainspan cortex scores for scatter plot
    """
    bs_scores_adult = bs_scores.loc['18 yrs':,:].groupby('structure_name').mean()

    both_scores = (pd.concat({'Brainspan': bs_scores_adult, 'AHBA': ahba_scores_mapped})
     .melt(ignore_index=False, var_name='G', value_name='score')
     .set_index('G', append=True)
     .unstack(0).droplevel(0, axis=1)
     .reset_index()
    )

    corrs = (both_scores
             .groupby('G').corr()
             .loc[(slice(None), 'AHBA'), 'Brainspan']
             .droplevel(1)
            )
    return both_scores, corrs


def make_brain_plots(version, version_to_bs_mapping, bs_scores):
    """
    Prepare brain scores for plotting comparison maps
    """
    ahba_scores_plot = (get_mapped_scores(version, version_to_bs_mapping, mean=False)
                        .loc[:,['G1','G2','G3']].reset_index()
                  )

    bs_scores_plot = (bs_scores
                      .loc['18 yrs':,:].groupby('structure_name').mean()
                      .apply(lambda x: (x-np.mean(x))/np.std(x))
                      .join(version_to_bs_mapping.set_index('structure_name')['label'])
                  )

    scores_plot = {
        'BrainSpan':bs_scores_plot,
        'AHBA':ahba_scores_plot
    }

    scores_plot = (
        pd.concat(scores_plot).reset_index(0).rename({'level_0':'version'},axis=1)
        .set_index(['version', 'label'])
        .reset_index()
    )
    
    return scores_plot


#########################################################



# Get AHBA and BS score comparison
def get_scores(version, bs_agg, hcp_bs, hcp_info, hcp_base):
    hcp_all = hcp_base.score_from(version).join(get_labels_hcp())
    
    hcp_filtered, hcp_mean = get_filtered_hcp(version, bs_agg, hcp_bs, hcp_info)
    
    bs_pcs_hcp = get_bs_pcs(version, bs_agg, hcp_bs)
    
    hcp_scores = (pd.concat({
        'AHBA_all': hcp_all,
        'AHBA_filtered': hcp_filtered.drop(['structure_name','cortex'],axis=1),
        'AHBA_mean': hcp_mean,
        'Brainspan': bs_pcs_hcp
    })
                  .reset_index(level=0).rename(columns={'level_0':'version'}).set_index(['version', 'label'])
                  .apply(lambda y: y.groupby('version').apply(lambda x: (x-np.mean(x))/np.std(x)), axis=0) # standardize all versions for plotting
                  .reset_index()
                 )
    
    cortex_scores = (hcp_scores
             .join(hcp_bs.set_index('region')['cortex'], on = 'label')
             .groupby(['version', 'cortex'])
             .mean()
             .set_axis([f'PC{i+1}' for i in range(5)], axis=1)
            )
    
    corrs = cortex_scores.loc['AHBA_mean'].corrwith(cortex_scores.loc['Brainspan'])

    cortex_scores = (cortex_scores
                     .loc[['AHBA_mean','Brainspan']]
                     .melt(ignore_index=False, var_name='PC')
                     .reset_index()
                     .pivot_table(index=['PC', 'cortex'],
                                  columns='version', values='value')
                     .loc[['PC1','PC2','PC3']].reset_index()
                    )
    
    return hcp_scores, cortex_scores, corrs

