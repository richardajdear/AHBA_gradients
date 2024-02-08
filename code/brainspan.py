# Helper functions for processing brainspan data

import numpy as np, pandas as pd
from processing import *
# import statsmodels.api as sm
from statsmodels.gam.api import GLMGam, BSplines
# import statsmodels.formula.api as smf
from itertools import combinations
from abagen.utils import efficient_corr


def get_brainspan_curves_by_gene(genes=None, return_points=False, alpha=0.5,
                                 save_path='../outputs/brainspan_curves.csv'):
    """
    1. Get BrainSpan data and match regions to HCP parcellation
    2. Convert age to continuous variable
    3. Fit and predict GAM curves for selected genes
    """
    if save_path is not None and genes is None:
        return pd.read_csv(save_path, index_col=0)
    
    # 1
    bs_exp, bs_col, bs_row = get_brainspan()
    hcp_bs_mapping = get_hcp_bs_mapping_v2()
    bs_clean = clean_brainspan(bs_exp, bs_col, bs_row, hcp_bs_mapping)
    # 2
    bs_continuous = make_continuous(bs_clean, norm_samples=True)
    # 3
    genes_to_fit = np.intersect1d(bs_clean.columns, genes)
    models = fit_gam_models(bs_continuous, genes_to_fit, alpha=alpha)
    curves, points = predict_gam_curves(models, bs_continuous, genes_to_fit)
    
    if save_path is not None:
        curves.to_csv(save_path)
        print(f"BrainSpan curves saved to {save_path}")

    if return_points:
        return curves, points
    else:
        return curves


def get_weight_quantiles(weights, q=10):
    quantiles = (weights
                 .melt(ignore_index=False, var_name='C',value_name='C_score')
                 .assign(C_quantile = lambda x: x.groupby('C')['C_score'].apply(lambda y: pd.qcut(y, q=q, labels=range(q))))
                 .rename_axis('gene').reset_index()
    )
    return quantiles

def get_quantile_curves(weights, curves, q=10):
    quantiles = get_weight_quantiles(weights)

    quantile_curves = (quantiles
     .join(curves.set_index('gene'), on='gene').dropna()
     .groupby(['C','C_quantile','age_log10'])
     .agg({'pred':'mean'})
     .reset_index()
    )
    return quantile_curves


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


def age_to_continuous(age_vector):
    """
    Parse BrainSpan age labels ('pcw','mos','yrs') to days post conception
    """
    age_vector = age_vector.astype('str')

    is_prenatal = age_vector.str.contains('pcw')
    pcw = age_vector[is_prenatal].str.replace(' pcw','').astype('int')

    is_months = age_vector.str.contains('mos')
    months = age_vector[is_months].str.replace(' mos','').astype('int')

    is_years = age_vector.str.contains('yrs')
    years = age_vector[is_years].str.replace(' yrs','').astype('int')

    age_vector[is_prenatal] = pcw*7 #-1*(40-pcw)/52
    age_vector[is_months] = 40*7+months*30 #months/12
    age_vector[is_years] = 40*7+years*365

    return age_vector.astype('float')


def make_continuous(bs_clean, log=True, norm_samples=True):
    """
    Convert BrainSpan age labels to continuous variable using age_to_continuous
    """
    if log:
        _bs = bs_clean.pipe(np.log10)
    else:
        _bs = bs_clean

    if norm_samples:
        _bs = _bs.apply(lambda x: x/np.mean(x), axis=1)

    return (_bs
        .rename_axis(None, axis=1)
        .reset_index()
        .assign(age = lambda x: age_to_continuous(x['age']))
        .assign(age_log10 = lambda x: np.log10(x['age']))
        .drop('structure_name', axis=1)
        .rename({'structure_acronym':'region'}, axis=1)
    )


def fit_gam_models(bs_continuous, genes_to_fit, age_var='age_log10', df=12, degree=3, alpha=1.0):
    """
    Fit GAM for each gene using continuous age
    """
    spline_x = bs_continuous[age_var]
    basis_splines = BSplines(spline_x, df=df, degree=degree)
    alpha=alpha

    models = {}
    for gene in genes_to_fit:
        # models[gene] = GLMGam(endog=bs_continuous[gene], exog=bs_continuous.loc[:, ['gender', 'region']],
        #                                 smoother=basis_splines, alpha=alpha).fit()
        formula = f'Q("{gene}") ~ {age_var} + gender + region'
        models[gene] = GLMGam.from_formula(formula=formula, data=bs_continuous, 
                                        smoother=basis_splines, alpha=alpha).fit()
        # models[gene] = smf.ols(formula=formula, data=bs_continuous).fit()
        
    return models


def predict_gam_curves(models, bs_continuous, genes_to_fit,
                       region='DFC', age_var='age_log10', n_preds=100):
    """
    Predict GAM curves for a single region and both genders
    """
    # Clean up input data
    df_data = (bs_continuous
                .drop('donor_id', axis=1)
                .loc[lambda x: x['region']==region, [age_var, 'gender'] + list(genes_to_fit)]
                .melt(id_vars=[age_var, 'gender'], var_name='gene', value_name='true')
                #  .assign(age = lambda x: (10**x[age_var]-40*7)/365)
    )

    # Vector of ages to predict at
    ages_to_predict = np.linspace(min(df_data[age_var]), max(df_data[age_var]), n_preds)
    # Dataframe to predict (repeat eachand F)
    df_preds = pd.DataFrame({
        age_var: np.repeat(ages_to_predict, 2),
        'gender': ['F','M'] * n_preds,
        'region': region
    })

    # Make predictions for each gene
    preds = {gene:models[gene].predict(df_preds, exog_smooth=df_preds[age_var]) for gene in genes_to_fit}
    # Combine genes into df
    df_preds = (df_preds
                .join(pd.concat(preds,axis=1))
                .melt(id_vars=[age_var, 'gender', 'region'], var_name='gene', value_name='pred')
                # Add gene predictions normalised to 75th quantile
                .assign(pred_q75=lambda x: x.groupby(['gene','gender','region'])
                        .apply(lambda y: y['pred']/np.quantile(y['pred'],.75)).reset_index([0,1,2], drop=True))
                # .assign(age = lambda x: (10**x['age_log10']-40*7)/365)
    )   
    return df_preds, df_data



def get_age_groups():
    """
    Define age groupings for use in Brainspan analysis
    """
    age_groups = {        
        # # Simple groupings ...
        '8 pcw': 'Pre-Birth',
        '9 pcw': 'Pre-Birth', #
        '12 pcw': 'Pre-Birth',
        '13 pcw': 'Pre-Birth',
        '16 pcw': 'Pre-Birth',
        '17 pcw': 'Pre-Birth',
        '19 pcw': 'Pre-Birth',
        '21 pcw': 'Pre-Birth',
        '24 pcw': 'Pre-Birth',
        '25 pcw': 'Pre-Birth', #
        '26 pcw': 'Pre-Birth', #
        '35 pcw': 'Pre-Birth', #
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
        '15 yrs': 'Birth-13 yrs', #
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

def count_samples(bs_col, bs_cortex_mapping):
    """
    Count samples from each Brainspan region and age
    """
    age_region_counts = (bs_col
                         .assign(cortex = lambda x: x['structure_name'].map(bs_cortex_mapping))
     .pivot_table(values='column_num', index='cortex', columns='age', aggfunc='count')
                        )
    return age_region_counts



def clean_brainspan(bs_exp, bs_col, bs_row, bs_mapping, log=False, norm_samples=False):
    """
    Clean up Brainspan data into dataframe
    Drop age groups with >=3 samples
    """
    mapped_regions = bs_mapping['structure_name'].dropna().unique()
    
    # Join region data
    # Filter for only mapped regions (i.e. no subcortex)
    bs_mapped = (pd.concat([
        bs_exp.T,
        bs_col.set_index('column_num')[['donor_id', 'gender', 'age', 'structure_name', 'structure_acronym']]
    ], axis=1)
     .loc[lambda x: np.isin(x['structure_name'], mapped_regions)] # filter for mapped regions
     .set_index(['donor_id', 'gender', 'age', 'structure_name', 'structure_acronym'])
     # .dropna(how='all')
    )
    
    # Join gene data and clean NAs and duplicates
    bs_clean = (bs_mapped
     .set_axis(bs_row['gene_symbol'], axis=1)
     .fillna(0)
     .loc[:, lambda x: (x != 0).all(axis=0)] # drop zero and na columns
     .loc[:, lambda x: ~x.columns.duplicated()] # some columns are duplicates
    )

    # Drop brains with < 4 samples
    bs_donor_counts = bs_clean.groupby(['donor_id']).size()
    bs_keep = bs_donor_counts > 3
    ix_keep = pd.IndexSlice[bs_keep[bs_keep].index.get_level_values('donor_id'), :]
    bs_clean = bs_clean.loc[ix_keep,:]

    if log:
        bs_clean = bs_clean.pipe(np.log10)

    if norm_samples:
        bs_clean = bs_clean.apply(lambda x: x/np.mean(x), axis=1)

    return bs_clean


def aggregate_brainspan_by_age(bs_clean, normalize=True):
    """
    Aggregate Brainspan donor brains by ages
    """
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
    bs_scores = (bs_agg.loc[:, gene_mask] @ version.weights.loc[gene_mask, :]).iloc[:, :3]
    
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
                       .melt(var_name='C', value_name='corr', ignore_index=False)
                       .reset_index())
    return bs_scores_corr


def combine_scores(bs_scores, ahba_scores_mapped, age_groups):
    """
    Combine AHBA and Brainspan cortex scores for scatter plot
    """
    bs_scores_adult = bs_scores.loc['18 yrs':,:].groupby('structure_name').mean()

    both_scores = (pd.concat({'Brainspan': bs_scores_adult, 'AHBA': ahba_scores_mapped})
     .melt(ignore_index=False, var_name='C', value_name='score')
     .set_index('C', append=True)
     .unstack(0).droplevel(0, axis=1)
     .reset_index()
    )

    corrs = (both_scores
             .groupby('C').corr()
             .loc[(slice(None), 'AHBA'), 'Brainspan']
             .droplevel(1)
            )
    return both_scores, corrs


def make_brain_plots(version, version_to_bs_mapping, bs_scores):
    """
    Prepare brain scores for plotting comparison maps
    """
    ahba_scores_plot = (get_mapped_scores(version, version_to_bs_mapping, mean=False)
                        .loc[:,['C1','C2','C3']].reset_index()
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



def get_stable_genes_brainspan(bs_clean, subjects_for_stability_test, 
                               percentile_threshold=0.5, return_stability=False):
    """
    Identify top X% most stable genes using differential stability
    Across selected subjects
    """
    # Get single subject expression for selected subjects
    subject_expression_dict = {
        s: bs_clean.loc[pd.IndexSlice[s,:,:,:,:]].reset_index(level=[0,1,2], drop=True) for s in subjects_for_stability_test
    }

    all_genes = np.array(bs_clean.columns)
    gene_corrs = np.zeros((len(all_genes), sum(range(len(subjects_for_stability_test)))))
    
    # Iterate over all pairs of subjects
    for n, (s1, s2) in enumerate(combinations(subjects_for_stability_test, 2)):
        matched_regions = np.intersect1d(subject_expression_dict[s1].index, subject_expression_dict[s2].index)
        gene_corrs[:, n] = efficient_corr(
                subject_expression_dict[s1].loc[matched_regions], 
                subject_expression_dict[s2].loc[matched_regions]
        )


    gene_stability = np.nanmean(gene_corrs, axis=1)
    threshold = np.nanpercentile(gene_stability, percentile_threshold*100)
    genes_to_keep = all_genes[gene_stability >= threshold]
    if return_stability:
        return gene_stability
    return genes_to_keep



## LEGACY


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




#########################################################



# # Get AHBA and BS score comparison
# def get_scores(version, bs_agg, hcp_bs, hcp_info, hcp_base):
#     hcp_all = hcp_base.score_from(version).join(get_labels_hcp())
    
#     hcp_filtered, hcp_mean = get_filtered_hcp(version, bs_agg, hcp_bs, hcp_info)
    
#     bs_pcs_hcp = get_bs_pcs(version, bs_agg, hcp_bs)
    
#     hcp_scores = (pd.concat({
#         'AHBA_all': hcp_all,
#         'AHBA_filtered': hcp_filtered.drop(['structure_name','cortex'],axis=1),
#         'AHBA_mean': hcp_mean,
#         'Brainspan': bs_pcs_hcp
#     })
#                   .reset_index(level=0).rename(columns={'level_0':'version'}).set_index(['version', 'label'])
#                   .apply(lambda y: y.groupby('version').apply(lambda x: (x-np.mean(x))/np.std(x)), axis=0) # standardize all versions for plotting
#                   .reset_index()
#                  )
    
#     cortex_scores = (hcp_scores
#              .join(hcp_bs.set_index('region')['cortex'], on = 'label')
#              .groupby(['version', 'cortex'])
#              .mean()
#              .set_axis([f'PC{i+1}' for i in range(5)], axis=1)
#             )
    
#     corrs = cortex_scores.loc['AHBA_mean'].corrwith(cortex_scores.loc['Brainspan'])

#     cortex_scores = (cortex_scores
#                      .loc[['AHBA_mean','Brainspan']]
#                      .melt(ignore_index=False, var_name='PC')
#                      .reset_index()
#                      .pivot_table(index=['PC', 'cortex'],
#                                   columns='version', values='value')
#                      .loc[['PC1','PC2','PC3']].reset_index()
#                     )
    
#     return hcp_scores, cortex_scores, corrs

