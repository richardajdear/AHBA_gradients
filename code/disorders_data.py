# Functions to read and clean data for disorders analysis

import numpy as np, pandas as pd
import mygene
import nibabel as nib
from neuromaps.nulls.spins import parcels_to_vertices, vertices_to_parcels
from processing import *
from scipy.stats import fisher_exact


replace_dict = {'01-Mar':'MARCH1', '02-Mar':'MARCH2', '03-Mar':'MARCH3', '04-Mar':'MARCH4', '05-Mar':'MARCH5', '06-Mar':'MARCH6', '07-Mar':'MARCH7', '08-Mar':'MARCH8', 
                '09-Mar':'MARCH9', '10-Mar':'MARCH10', '11-Mar':'MARCH11',
                '01-Sep':'SEPT1', '02-Sep':'SEPT2', '03-Sep':'SEPT3', '04-Sep':'SEPT4', '05-Sep':'SEPT5', '06-Sep':'SEPT6', '07-Sep':'SEPT7', '08-Sep':'SEPT8',
                '09-Sep':'SEPT9', '10-Sep':'SEPT10', '11-Sep':'SEPT11', '12-Sep':'SEPT12', '13-Sep':'SEPT13', '14-Sep':'SEPT14', '15-Sep':'SEPT15', 
                '01-Dec':'DECR1', '02-Dec':'DECR2'}


def ensembl_id_to_gene_symbol(ensembl_ids):
    mg = mygene.MyGeneInfo()
    ensembl_ids = [ens.split('.')[0] for ens in ensembl_ids]
    matches = mg.querymany(ensembl_ids, scopes='ensembl.gene', as_dataframe=True)['symbol']
    return(matches)


def get_disorder_maps(data_dir="../data/lifespan_dx_DKatlas.csv"):
    maps = (
        pd.read_csv(data_dir, index_col=0)
        # .apply(lambda x: (x-np.mean(x))/np.std(x))
        .sort_index(axis=1)
        .rename_axis('label')#.reset_index()
    )
    maps = maps.set_index(maps.index.str.replace(" ", ""))
 
    maps = maps.set_index('lh_' + maps.index)
    return maps


def get_gwas_combined():
    # SCZ data from https://figshare.com/articles/dataset/scz2022/19426775?file=35775617
    trubetskoy = (pd.read_csv(f"../data/gwas/trubetskoy2022_extended.csv")
                #   .loc[lambda x: x["Extended.GWAS"]=='YES', 'Symbol.ID']
                #   .rename('gene')
                  .loc[lambda x: x["Extended.GWAS"]=='YES', 'Ensembl.ID']
                  .pipe(ensembl_id_to_gene_symbol)
                  .rename('gene')
                  .reset_index(drop=True)
    )

    # ASD data from https://www.nature.com/articles/s41398-020-00953-9#MOESM1
    matoba = (pd.read_csv(f"../data/gwas/matoba2020_tableS7.csv")
                  .loc[lambda x: x['FDR']<=0.05, :]
                #   .loc[:, 'hgnc_symbol'].rename('gene')
                  .loc[:, 'GENE']
                  .pipe(ensembl_id_to_gene_symbol)
                  .rename('gene')
                  .reset_index(drop=True)
    )

    # MDD data from https://www.nature.com/articles/s41593-018-0326-7#MOESM11
    howard = (pd.read_csv(f"../data/gwas/howard2019_tableS9.csv", header=1)
                  .loc[:, 'Gene Name'].rename('gene')
    )

    df = (pd.concat({
        'ASD': matoba,
        'MDD': howard,
        'SCZ': trubetskoy,
    })
        .reset_index(0)
        .rename({'level_0':'label'}, axis=1)
        .assign(gene = lambda x: x['gene'].str.replace('\\..*','', regex=True)) #drop variants
        .replace({'gene': replace_dict})
        .drop_duplicates()
        .dropna()
    )
    return df




def get_deg_combined(scz_only=False, updown=False, logFC = False):
    deg_dict = {
        'ASD--Gandal 2022': pd.read_csv("../data/deg/gandal2022_tableS3.csv").loc[lambda x: x['WholeCortex_ASD_FDR']<=0.05, ['external_gene_name', 'WholeCortex_ASD_logFC']],
        # 'ASD--Ramaswami 2020': pd.read_csv("../data/deg/ramaswami2020_tableS2.csv", header=1).loc[lambda x: x['p.condition.fdr']<=0.05, 'Gene'].pipe(ensembl_id_to_gene_symbol),
        'ASD--Gandal 2018': pd.read_csv("../data/deg/gandal_genes_rnaseq.csv").loc[lambda x: x['ASD.fdr']<=0.05, ['gene_name', 'ASD.log2FC']],
        'ASD--Parikshak 2016': pd.read_csv("../data/deg/parikshak2016_tableS2.csv", header=1, dtype={'Chr':'object', 'HGNC Symbol':'object'}).loc[lambda x: x['FDR-adjusted P value, ASD vs CTL']<=0.05, ['HGNC Symbol', 'log2(FC) ASD vs CTL']],
        'MDD--Jaffe 2022': pd.read_csv("../data/deg/jaffe2022_dataS1.csv", header=0, usecols=['Symbol','Cortex_logFC_MDD','Cortex_adjPVal_MDD']).loc[lambda x: x['Cortex_adjPVal_MDD']<=0.05, ['Symbol','Cortex_logFC_MDD']],
        # 'MDD--Girgenti 2021': pd.read_csv("../data/deg/girgenti2021_tableS9.csv", header=0, usecols=['Genename','MDD.dlPFC.log2FoldChange','MDD.dlPFC.padj']).loc[lambda x: x['MDD.dlPFC.padj']<=0.05, ['Genename','MDD.dlPFC.log2FoldChange']],
        'SCZ--Gandal 2018': pd.read_csv("../data/deg/gandal_genes_rnaseq.csv").loc[lambda x: x['SCZ.fdr']<=0.05, ['gene_name', 'SCZ.log2FC']],
        'SCZ--Fromer 2016': pd.read_csv("../data/deg/fromer2016_tableS3.csv", header=1).loc[lambda x: x['FDR estimate']<=0.05, ['Gene Symbol','logFC']],
        'SCZ--Collado-Torres 2019': pd.read_csv("../data/deg/colladotorres2019_tableS11.csv").query("region=='DLPFC'").loc[lambda x: x['adj.P.Val']<=0.05, ['Symbol', 'logFC']],
        'SCZ--Jaffe 2018': pd.read_csv("../data/deg/jaffe2018_tableS9.csv", header=1).loc[lambda x: x['fdr_qsva']<=0.05, ['Symbol', 'log2FC_qsva']],
    }
    deg_dict = {name:deg.set_axis(['gene', 'log2FC'], axis=1) for name, deg in deg_dict.items()}

    deg_genes = (pd.concat(deg_dict)
                 .reset_index(0).rename({'level_0':'label'}, axis=1)
                 .replace({'gene': replace_dict})
    )

    if updown:
        deg_genes = deg_genes.assign(updown = lambda x: np.where(x['log2FC']>0,'up','down'))

    if not logFC:
        # drop logFC to cleanly remove duplicates
        deg_genes = deg_genes.drop('log2FC', axis=1)

    deg_genes = (deg_genes
                 .drop_duplicates()
                 .dropna()
    )

    if scz_only:
        deg_genes = (deg_genes
                     .loc[lambda x: x['label'].str.contains('SCZ'),:]
                     .assign(label = lambda x: x['label'].str.replace("SCZ--",""))
        )

    deg_genes = deg_genes.assign(label = lambda x: pd.Categorical(x['label'], 
                                                    ordered=True, categories=x['label'].unique()))
    return deg_genes


def get_deg_consensus(updown=False, logFC=False):
    deg_all = get_deg_combined(updown=updown, logFC=logFC)

    # Optionally take only consensus with directional agreement
    if updown:
        grouping = ['label','gene','updown']
    else:
        grouping = ['label','gene']

    if logFC:
        agg_dict = {'study':'nunique', 'log2FC': 'mean'}
    else:
        agg_dict = {'study':'nunique'}

    deg_consensus = (deg_all
                    .assign(
                        study = lambda x: x['label'].str.split('--', expand=True)[1],
                        label = lambda x: x['label'].str.split('--', expand=True)[0]
                        )
                    .groupby(grouping, as_index=False).agg(agg_dict)
                    .loc[lambda x: (x['study'] >= 2) | (x['label']=='MDD')] # either 2+ studies, or MDD
                    )
    
    return deg_consensus



def get_scz_gyral_sulcal(
            data_path = '../data/lh.GyralSulcalDifferences.mgh', 
            hcp_img = '../data/parcellations/lh.HCPMMP1.annot',
            name = 'SCZ supragranular'
        ):
    scz_diff_img = nib.load(data_path)
    hcp_img = annot_to_gifti(hcp_img)
    
    scz_diff_hcp = vertices_to_parcels(scz_diff_img.get_fdata(), hcp_img)

    scz_diff = pd.Series(scz_diff_hcp,
                         index = get_labels_hcp()[:180],
                         name = name
                         )

    scz_diff = (scz_diff
                .to_frame()
                .apply(lambda x: (x-np.mean(x))/np.std(x))
                .apply(lambda x: -x) # invert to focus on thinning
                )

    return scz_diff


def get_pergola_consensus():
    return [
        'CHRNA4',
        'DENND1A',
        'ADCY9',
        'ADCY5',
        'MKL1',
        'DLGAP2',
        'RTL10',
        'HECW1',
        'RAB3A',
        'ATP2A2',
        'UBE2O',
        'KCNQ3',
        'NOS1AP',
        'KIF3C',
        'DPP6',
        'SRRM4',
        'KSR2',
        'SLITRK1',
        'OPCML',
        'PI4KA',
        'KIAA1549L',
        'UNC79',
        'CLSTN3',
        'UNC80',
        'TENM2',
        'FAM135B',
        'PTPRN2',
        'RGAG',
    ]

### LEGACY


# def get_disgenet_genes():
#     disgenet_genes = (pd.read_csv("../data/disgenet_genes.csv")
#         .assign(score_pct = lambda x: x.groupby('Disease')['Score_gda'].rank(pct=True))
#         # .loc[lambda x: x['score_pct'] >= 0.75]
#         .rename({'Disease':'label','Gene':'gene'},axis=1)
#         .replace({
#             'Autism Spectrum Disorders':'ASD',
#             'Major Depressive Disorder':'MDD',
#             'Schizophrenia':'SCZ'
#         })
#     )
#     return disgenet_genes 


# def get_gwas(filename='trubetskoy2022'):
#     df = pd.read_csv(f"../data/gwas/{filename}.csv")
#     df = (df.melt(var_name='label', value_name='gene').dropna()
#           .assign(gene = lambda x: x['gene'].str.replace('\\..*','', regex=True)) #drop variants
#           .drop_duplicates()
#           .replace({'gene': replace_dict})
#     )
#     return df




# def get_asd_deg():
#     asd_deg = (pd.concat({
#         'gandal2022': pd.read_csv("../data/deg/gandal2022_asd.csv").dropna(),
#         'parikshak2016': pd.read_csv("../data/deg/parikshak2016_genes.csv").dropna(),
#     }).reset_index(0).rename({'level_0':'label'},axis=1)
#         .loc[lambda x: x['FDR']<0.05]
#         .drop('FDR',axis=1)
#         .replace({'gene': replace_dict})
#     )

#     gandal2018 = (
#         get_gandal_genes(which='rnaseq', disorders=['ASD'], sig_level=0.05)
#         .reset_index().query('sig').drop('sig', axis=1)
#         .assign(label='gandal2018')
#     )

#     return pd.concat([asd_deg, gandal2018])




# def get_gandal_genes(which='microarray', disorders = ['ASD', 'SCZ', 'MDD'], sig_level=.05):

#     if which == 'microarray':
#         gandal_genes = (pd.read_csv("../data/deg/gandal_genes_microarray.csv")
#                         .rename({'external_gene_id':'gene'},axis=1)
#                         .replace({'gene': replace_dict})
#                         .rename({f'{d}.beta_log2FC':f'{d}.log2FC' for d in disorders}, axis=1) 
#                         .set_index('gene')
#                         .assign(rank = lambda x: x.groupby('gene')['start_position'].rank(method='first'))
#                         .loc[lambda x: x['rank']==1]
#                         # .loc[lambda x: x.groupby('gene')['rank'].max()>1]
#                         # .sort_index()
#                        )
#     elif which == 'rnaseq':
#         gandal_genes = (pd.read_csv("../data/deg/gandal_genes_rnaseq.csv")
#                         .rename({'gene_name':'gene'},axis=1) 
#                         .replace({'gene':replace_dict})
#                         .rename({f'{d}.fdr':f'{d}.FDR' for d in disorders}, axis=1)                         
#                         .set_index('gene')
#                         .assign(rank = lambda x: x.groupby('gene')['transcript_length'].rank(method='first'))
#                         .loc[lambda x: x['rank']==1]
#                         # .loc[lambda x: x.groupby('gene')['rank'].max()>1]
#                         # .sort_index()
#                        )
        
#     for d in disorders:
#         gandal_genes[f'{d}.sig'] = gandal_genes[f'{d}.FDR'] < sig_level

#     gandal_sig = (gandal_genes
#             .loc[:,[f'{d}.sig' for d in disorders]]
#             .melt(ignore_index=False, var_name='label', value_name='sig')
#             .assign(label = lambda x: x['label'].str.replace('.sig','', regex=False))
#             .set_index('label', append=True)
#             )

#     gandal_log2FC = (gandal_genes
#             .loc[:,[f'{d}.log2FC' for d in disorders]]
#             .melt(ignore_index=False, var_name='label', value_name='log2FC')
#             .assign(label = lambda x: x['label'].str.replace('.log2FC','', regex=False))
#             .set_index('label', append=True)
#             .join(gandal_sig)
#         )

#     return gandal_log2FC


# def get_gene_weights_with_labels(weights, gandal_genes, 
#                                  disorders = ['ASD','SCZ','MDD']):
#     gandal_fdr = (gandal_genes
#                   .loc[:,[f'{d}.FDR' for d in disorders]]
#                   .melt(ignore_index=False, var_name='disorder', value_name='FDR')
#                   # .loc[lambda x: x['FDR']]
#                   .assign(disorder = lambda x: 
#                           x['disorder'].str.replace('.FDR','', regex=False))
#                   .set_index('disorder', append=True)
#                  )

#     gandal_log2FC = (gandal_genes
#                   .loc[:,[f'{d}.log2FC' for d in disorders]]
#                   .melt(ignore_index=False, var_name='disorder', value_name='log2FC')
#                   .assign(disorder = lambda x: 
#                           x['disorder'].str.replace('.log2FC','', regex=False))
#                   .set_index('disorder', append=True)
#                   .join(gandal_fdr)
#                   .reset_index('disorder')
#                   .loc[lambda x: x['FDR'] < 0.05]
#                  )

#     weights_labels = (weights
#                       .join(gandal_log2FC)     
#                       .rename({'disorder':'label'}, axis=1)                   
#                       .fillna({'label':'none'})
#                       .assign(cluster = lambda x: 
#                               np.select([x['log2FC']>0,x['log2FC']<0],
#                                         ['upregulated','downregulated'], 
#                                         default=None)
#                              )
#                      )
#     return weights_labels


# def get_rare_genes():
#     rare_genes = pd.read_csv("../data/rare_genes_konrad.csv").melt(var_name='label',value_name='gene').dropna()
#     return rare_genes




# def get_gene_corr(weights, null_weights, gandal_genes, sig_thresh=1, adjust='fdr_bh', disorders = ['ASD', 'MDD', 'SCZ']):
#     """
#     Correlate gene weights with log2FC
#     """
#     true_corrs = {}
#     null_corrs = {}

#     for d in disorders:
#         # Skip if disorder not present
#         if f'{d}.log2FC' not in gandal_genes.columns:
#             continue

#         # Don't overwrite the input data in the loop!
#         genes_log2FC = gandal_genes.copy()            

#         # Take only significant genes
#         genes_log2FC = genes_log2FC.loc[lambda x: x[f'{d}.FDR']<sig_thresh, f'{d}.log2FC']
        
#         # Find matching genes and filter differential expression
#         genes_matched = set(weights.index).intersection(genes_log2FC.index)
#         genes_log2FC = genes_log2FC.loc[genes_matched].values
#         # Also filter true and null weights
#         # np.searchsorted to find indices of matched genes in same order as genes_matched
#         genes_idxs = np.searchsorted(weights.index, list(genes_matched))
#         true = weights.values[genes_idxs, :]
#         nulls = null_weights[genes_idxs, :, :]
        
#         # For each PC, get the correlation of the true weights and all the nulls
#         true_corr_disorder = np.zeros(3)
#         null_corr_disorder = np.zeros((null_weights.shape[2], 3))
#         for i in range(3):
#             true_corr_disorder[i] = np_pearson_corr(genes_log2FC, true[:,i]).squeeze()
#             null_corr_disorder[:, i] = np_pearson_corr(genes_log2FC, nulls[:,i,:]).squeeze()
        
#         true_corrs[d] = pd.DataFrame(true_corr_disorder).T
#         null_corrs[d] = pd.DataFrame(null_corr_disorder) 

#     true_corrs = pd.concat(true_corrs).set_axis(['G1','G2','G3'], axis=1).droplevel(1)
#     null_corrs = pd.concat(null_corrs).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)

#     null_p = compute_null_p(true_corrs, null_corrs, adjust=adjust)
#     return null_p



# def get_gene_corr2(weights, null_weights, genes_log2FC, group='disorder', sig_thresh=1, adjust='fdr_bh'):
#     """
#     Correlate gene weights with log2FC
#     """
#     true_corrs = {}
#     null_corrs = {}

#     for g in genes_log2FC[group].unique():
#         # Don't overwrite the input data in the loop!
#         _genes_log2FC = genes_log2FC.copy().loc[lambda x: x[group]==g]

#         # Take only significant genes
#         _genes_log2FC = _genes_log2FC.loc[lambda x: x['FDR']<sig_thresh, 'log2FC']
        
#         # Merge duplicates
#         _genes_log2FC = _genes_log2FC.groupby('gene').mean()

#         # Find matching genes and filter differential expression
#         genes_matched = list(set(weights.index).intersection(_genes_log2FC.index))
#         _genes_log2FC = _genes_log2FC.loc[genes_matched].values
#         # Also filter true and null weights
#         # np.searchsorted to find indices of matched genes in same order as genes_matched
#         genes_idxs = np.searchsorted(weights.index, list(genes_matched))
#         true = weights.values[genes_idxs, :]
#         nulls = null_weights[genes_idxs, :, :]
        
#         # For each PC, get the correlation of the true weights and all the nulls
#         true_corr_group = np.zeros(3)
#         null_corr_group = np.zeros((null_weights.shape[2], 3))
#         for i in range(3):
#             true_corr_group[i] = np_pearson_corr(_genes_log2FC, true[:,i]).squeeze()
#             null_corr_group[:, i] = np_pearson_corr(_genes_log2FC, nulls[:,i,:]).squeeze()
        
#         true_corrs[g] = pd.DataFrame(true_corr_group).T
#         null_corrs[g] = pd.DataFrame(null_corr_group) 

#     true_corrs = pd.concat(true_corrs).set_axis(['G1','G2','G3'], axis=1).droplevel(1)
#     null_corrs = pd.concat(null_corrs).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)

#     # return true_corrs, null_corrs
#     null_p = compute_null_p(true_corrs, null_corrs, adjust=adjust)
#     return null_p




# def get_gene_map_corr(version, null_scores, gandal_genes, posneg, method='mean', sig_thresh=1, 
#                       disorders = ['ASD','MDD','SCZ'], return_maps=False):
    
#     disorders = [d for d in disorders if f'{d}.log2FC' in gandal_genes.columns]

#     gene_log2FC = gandal_genes.copy()
    
#     # Set non sig changes to 0
#     for d in disorders:
#         gene_log2FC[f'{d}.log2FC'] = np.where(gene_log2FC[f'{d}.FDR'] >= sig_thresh, 0, gene_log2FC[f'{d}.log2FC'])

#     gene_log2FC = (gene_log2FC
#                      .loc[:, [f'{d}.log2FC' for d in disorders]]
#                      .rename({f'{d}.log2FC':f'{d}' for d in disorders}, axis=1)
#                      .melt(ignore_index=False, var_name='label', value_name='weight')
#                      .reset_index()
#                     )

#     if posneg == 'pos':
#         gene_log2FC['weight'] = gene_log2FC['weight'].clip(lower=0)
#     elif posneg == 'neg':
#         gene_log2FC['weight'] = gene_log2FC['weight'].clip(upper=0)*-1
#     elif posneg == 'abs':
#         gene_log2FC['weight'] = gene_log2FC['weight'].abs()

#     gene_maps = make_gene_maps(version, gene_log2FC, normalize='std', method=method)
    
#     if return_maps:
#         return gene_maps
    
#     true_corrs = get_corrs(version.clean_scores(), gene_maps)
#     null_corrs = corr_nulls_from_grads(null_scores, version.clean_scores(), gene_maps)
#     null_corrs_p = get_null_p(true_corrs, null_corrs, adjust='fdr_bh').rename({'map':'label'},axis=1)
    
#     return null_corrs_p


# def get_gene_sig(weights, null_weights, gandal_genes, posneg=None, posneg_weights=None, disorders = ['ASD', 'MDD', 'SCZ'], combine=None):
#     genes_sig_dict = {}
    
#     for d in disorders:
#         # Skip if disorder not present
#         if f'{d}.log2FC' not in gandal_genes.columns:
#             continue
        
#         # Don't overwrite the input data in the loop!
#         gene_sig = gandal_genes.copy()
        
#         # Filter for up or down regulated
#         if posneg == 'pos':
#             gene_sig = gene_sig.loc[lambda x: x[f'{d}.log2FC'] > 0]
#         elif posneg == 'neg':
#             gene_sig = gene_sig.loc[lambda x: x[f'{d}.log2FC'] < 0]
#         elif posneg == 'abs':
#             pass
        
#         # Filter for significant only
#         genes_sig = gene_sig.loc[lambda x: x[f'{d}.FDR'] < .05]
#         genes_sig_dict[d] = pd.Series(genes_sig.index)

#     # Combine across disorders
#     genes_sig = pd.concat(genes_sig_dict).reset_index(0).rename({'gene_name':'gene','level_0':'label'},axis=1)
    
#     # Take superset if desired
#     if combine=='union':
#         genes_sig['label'] = '-'.join(genes_sig['label'].unique())
#         genes_sig = genes_sig.drop_duplicates()
#     if combine=='ASD-SCZ intersection':
#         ASD_genes = genes_sig.loc[lambda x: x['label'] == 'ASD', 'gene']
#         SCZ_genes = genes_sig.loc[lambda x: x['label'] == 'SCZ', 'gene']
#         intersection = set(ASD_genes).intersection(SCZ_genes)
#         genes_sig = genes_sig.set_index('gene').loc[intersection, :].assign(label = 'ASD-SCZ\nintersection').reset_index()
    
#     true_mean, null_mean = compute_enrichments(weights, null_weights, genes_sig, posneg=posneg_weights)
#     null_p = compute_null_p(true_mean, null_mean, adjust='fdr_bh')

#     return null_p


# def get_gene_dot(weights, null_weights, gandal_genes, posneg=None, posneg_weights=None, sig_thresh=1, disorders = ['ASD', 'MDD', 'SCZ']):
#     """
#     x
#     """
#     weights = weights.copy()
#     null_weights = null_weights.copy()
#     if posneg_weights =='abs':
#         weights = np.abs(weights)
#         null_weights = np.abs(null_weights)
#     elif posneg_weights=='pos':
#         weights = pd.DataFrame(np.where(weights<0, 0, weights), index=weights.index)
#         null_weights = np.where(null_weights<0, 0, null_weights)
#     elif posneg_weights=='neg':
#         weights = pd.DataFrame(np.where(weights>0, 0, weights), index=weights.index)
#         null_weights = np.where(null_weights>0, 0, null_weights)    
    
#     true_dot = {}
#     null_dot = {}
#     for d in disorders:
#         # Skip if disorder not present
#         if f'{d}.log2FC' not in gandal_genes.columns:
#             continue          
            
#         # Don't overwrite the input data in the loop!
#         genes_log2FC = gandal_genes.copy()           
            
#         # Take only significant genes
#         genes_log2FC = genes_log2FC.loc[lambda x: x[f'{d}.FDR']<sig_thresh, f'{d}.log2FC']

#         # Find matching genes and filter differential expression
#         genes_matched = set(weights.index).intersection(genes_log2FC.index)
#         genes_log2FC = genes_log2FC.loc[genes_matched].values
#         # Also filter true and null weights
#         # np.searchsorted to find indices of matched genes in same order as genes_matched
#         genes_idxs = np.searchsorted(weights.index, list(genes_matched))
#         true = weights.values[genes_idxs, :]
#         nulls = null_weights[genes_idxs, :, :]
        
#         # Clip or take absolute
#         if posneg == 'pos':
#             genes_log2FC = genes_log2FC.clip(min=0)
#         elif posneg == 'neg':
#             genes_log2FC = genes_log2FC.clip(max=0)*-1
#         elif posneg == 'abs':
#             genes_log2FC = genes_log2FC.abs()

#         # Swap nulls axes for dot product
#         nulls = np.swapaxes(nulls[:, :, :], 0, 2)
                
#         # Compute dot product
#         true_dot[d] = pd.Series(true.T @ genes_log2FC)
#         null_dot[d] = pd.DataFrame(np.dot(nulls, genes_log2FC))

#     true_dot = pd.concat(true_dot).unstack(1).set_axis(['G1','G2','G3'], axis=1)
#     null_dot = pd.concat(null_dot).set_axis(['G1','G2','G3'], axis=1).reset_index(level=0).rename({'level_0':'label'}, axis=1)
    
#     null_p = compute_null_p(true_dot, null_dot, adjust='fdr_bh')
    
#     return null_p


