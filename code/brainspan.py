# Helper functions for processing brainspan data

import pandas as pd
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


# Define
def get_age_groups():
    age_groups = {
        '8 pcw': '8-12 pcw',
        '9 pcw': '8-12 pcw',
        '12 pcw': '8-12 pcw',
        '13 pcw': '13-18 pcw',
        '16 pcw': '13-18 pcw',
        '17 pcw': '13-18 pcw',
        '19 pcw': '19-24 pcw',
        '21 pcw': '19-24 pcw',
        '24 pcw': '19-24 pcw',
        '25 pcw': '25-38 pcw',
        '26 pcw': '25-38 pcw',
        '35 pcw': '25-38 pcw',
        '37 pcw': '25-38 pcw',
        '4 mos': 'Birth-5 months',
        '10 mos': '6-18 months',
        '1 yrs': '6-18 months',
        '2 yrs': '19 months-5yrs',
        '3 yrs': '19 months-5yrs',
        '4 yrs': '19 months-5yrs',
        '8 yrs': '6-11 yrs',
        '11 yrs': '6-11 yrs',
        '13 yrs': '12-19 yrs',
        '15 yrs': '12-19 yrs',
        '18 yrs': '12-19 yrs',
        '19 yrs': '12-19 yrs',
        '21 yrs': '20-40 yrs',
        '23 yrs': '20-40 yrs',
        '30 yrs': '20-40 yrs',
        '36 yrs': '20-40 yrs',
        '37 yrs': '20-40 yrs',
        '40 yrs': '20-40 yrs',   
    }
    return age_groups


# Filter HCP version by brainspan
def get_filtered_hcp(version, bs_agg, hcp_bs, hcp_info):
    hcp_filtered = (version.scores
     .join(get_labels_hcp())
     .rename_axis('id')
     .join(hcp_bs.set_index('region'), on='label')
     .dropna(axis=0)
    )

    hcp_mean = (hcp_filtered
     .groupby('cortex')
     .mean()
     .join(hcp_info.set_index('cortex')[['region']])
     .rename(columns={'region':'label'})
    )
    
    return hcp_filtered, hcp_mean

# Get BS PCs in HCP
def get_bs_pcs(version, bs_agg, hcp_bs):
    gene_mask = version.coefs.columns.intersection(bs_agg.columns)

    # Compute PCs
    bs_pcs = (bs_agg.loc[:, gene_mask] @ version.coefs.loc[:, gene_mask].T).iloc[:, :5]

    # Project into HCP regions
    bs_pcs_hcp = (bs_pcs
                  .reset_index(level=0, drop=True)
                  .join(hcp_bs.set_index('structure_name')['region'])
                  .rename(columns={'region':'label'})
                 ) 
    
    return bs_pcs_hcp

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

