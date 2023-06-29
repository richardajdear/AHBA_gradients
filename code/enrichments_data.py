# Functions to read in external data for enrichments analysis
import numpy as np, pandas as pd

# Dictionary of gene labels that get overwritten to dates in excel...
replace_dict = {'01-Mar':'MARCH1', '02-Mar':'MARCH2', '03-Mar':'MARCH3', '04-Mar':'MARCH4', '05-Mar':'MARCH5', '06-Mar':'MARCH6', '07-Mar':'MARCH7', '08-Mar':'MARCH8', 
                '09-Mar':'MARCH9', '10-Mar':'MARCH10', '11-Mar':'MARCH11',
                '01-Sep':'SEPT1', '02-Sep':'SEPT2', '03-Sep':'SEPT3', '04-Sep':'SEPT4', '05-Sep':'SEPT5', '06-Sep':'SEPT6', '07-Sep':'SEPT7', '08-Sep':'SEPT8',
                '09-Sep':'SEPT9', '10-Sep':'SEPT10', '11-Sep':'SEPT11', '12-Sep':'SEPT12', '13-Sep':'SEPT13', '14-Sep':'SEPT14', '15-Sep':'SEPT15', 
                '01-Dec':'DECR1', '02-Dec':'DECR2'}


### GO enrichments from STRING db

def clean_go_enrichment(file, direction=None, clean_terms=True, FDR_filter=None):
    """
    Clean STRING enrichment file
    Combine 'Antigen processing' enrichments into one
    """
    enrichment = pd.read_csv(file, delimiter='\t')
    
    if clean_terms:
        # Combine antigen terms
        enrichment.replace({'term description': "Antigen processing.*"}, 'Antigen processing', regex=True, inplace=True)
        # enrichment.drop_duplicates('term description', inplace=True)
        # Shorten other names
        replacements = {
            'Negative regulation of cytokine production involved in immune response': 'Regulation of cytokine production',
            'Negative regulation of natural killer cell mediated cytotoxicity': 'Regulation of killer cell mediated cytotoxicity',
            'Protection from natural killer cell mediated cytotoxicity': 'Protection from killer cell mediated cytotoxicity',
            'protein targeting to membrane':'protein targeting',
            'Positive regulation of extrinsic apoptotic signaling pathway': 'Regulation of apoptotic signaling pathway',
            'multicellular organism': '',
            'Establishment of protein localization to': 'Protein localization to',
            'Nuclear-transcribed mrna catabolic process,.*':'mRNA catabolic process',
            'Energy derivation by oxidation of organic compounds':'Energy derivation by oxidation',
            'Negative regulation of gene expression, epigenetic':'Regulation of gene expression, epigenetic',
            'Cellular process involved in reproduction in':'Cellular process involved in reproduction',
            'Regulation of dendritic spine development':'Dendritic spine development',
            'Regulation of dendrite development':'Dendrite development',
            'Serotonin receptor signaling pathway':'Serotonin receptor signaling',
            'SRP-dependent cotranslational protein targeting': 'SRP-dependent protein targeting',
            'Cellular response to catecholamine stimulus': 'Response to catecholamine',
            'Positive regulation of neurotransmitter secretion': 'Regulation of neurotransmitter secretion',
            'Cellular response to dopamine': 'Response to dopamine',
            'Cellular response to jasmonic acid stimulus': 'Response to jasmonic acid',
            'Negative regulation of chromatin silencing': 'Regulation of chromatin silencing',
            'Regulation of mrna splicing, via spliceosome': 'Regulation of mRNA splicing',
            'Negative regulation of mitochondrial fusion': 'Regulation of mitochondrial fusion',
            'Generation of precursor metabolites and energy': 'Generation of metabolites and energy',
            'Regulation of short-term neuronal synaptic plasticity': 'Regulation of short-term plasticity',
            'Regulation of macrophage migration.*': 'Regulation of macrophage migration',
            'Regulation of gene expression, epigenetic': 'Regulation of epigenetic expression'

        }
        enrichment.replace(replacements, inplace=True, regex=True)
        
        enrichment.drop_duplicates('term description', inplace=True)
    
    if direction is not None:
        enrichment = enrichment.loc[lambda x: x['direction'] == direction]
    else:
        enrichment = enrichment.loc[lambda x: x['direction'] != 'both ends']
    
    if FDR_filter is not None:
        enrichment = enrichment.loc[lambda x: x['false discovery rate'] <= FDR_filter]
    
    enrichment = (enrichment
                  .assign(FDR = lambda x: x['false discovery rate'],
                          neglogFDR = lambda x: -np.log10(x['false discovery rate']),
                          enrichment = lambda x: x['enrichment score'],
                          n_genes = lambda x: x['genes mapped'],
                          description = lambda x: x['term description']
                         )
                  .loc[:, ['description', 'n_genes', 'direction', 'enrichment', 'FDR', 'neglogFDR']]
                 )
    return enrichment


def combine_go_enrichments(version_, type_, dir_="../outputs/string_data/", 
                        include_c1=False, directed=False, top_n=None, **kwargs):
    """
    Combine enrichments from STRING
    """
    directions = {'C1':'top', 'C2':'bottom', 'C3':'bottom'}
    if not directed:
        directions = {k:None for k in directions.keys()}
    
    d = {
        'C2': clean_go_enrichment(dir_ + version_ + '/g2' + '.enrichment.' + type_ + '.tsv', direction = directions['C2'], **kwargs),
        'C3': clean_go_enrichment(dir_ + version_ + '/g3' + '.enrichment.' + type_ + '.tsv', direction = directions['C3'], **kwargs)
    }
    
    if include_c1:
        d.update({
            'C1': clean_go_enrichment(dir_ + version_ + '/g1' + '.enrichment.' + type_ + '.tsv', direction = directions['C1'], **kwargs)
        })
    
    df = (pd.concat(d)
          .reset_index(0).rename({'level_0':'C'}, axis=1)
          .assign(FDR_rank = lambda x: x.groupby(['C','direction'])['neglogFDR'].rank(method='first', ascending=False))
          .assign(enrichment_rank = lambda x: x.groupby(['C','direction'])['enrichment'].rank(method='first', ascending=False))
         )
    
    if directed:
        df['direction'] = 'bottom'
    
    if top_n is not None:
        df = (df
              .assign(rank = lambda x: x.groupby(['C','direction'])['neglogFDR']
                                        .rank(method='first', ascending=False))
              .loc[lambda x: x['rank'] <= top_n]
              .set_index(['C','direction'])
              .assign(maxrank = lambda x: x.groupby(['C','direction'])['rank'].max())
              .assign(rank = lambda x: x['maxrank']-x['rank'] + 1)
              .reset_index()
            #   .assign(rank = lambda x: top_n+1-x['rank'])
            )
    
    return df


def get_cell_genes(which='jakob', include=None, subtype=False, combine_layers=False, combine_ex_in=False, add_synapses=True):
    """
    Read cell genes table
    """

    if which == 'wen':
        path = '../data/wen_cell_genes.csv'
        cell_genes = (
            pd.read_csv(path)
            .set_axis(['gene','label'],axis=1)
            .loc[lambda x: np.isin(x['label'], ['purple','brown','blue'])]
            .replace({'label':{
                'purple':'Neurons', 'brown':'Synapses', 'blue':'Glia'
            }})
            .set_index('label').loc[['Neurons','Synapses','Glia'],:].reset_index()
        )
    elif which == 'zeng':
        path = '../data/zeng_layers.csv'
        cell_genes = pd.read_csv(path).set_axis(['label','gene'], axis=1)
    elif which == 'jakob':
        path="../data/jakob_cell_genes.csv"
        cell_genes = pd.read_csv(path)

        if include == 'only_lake':
            cell_genes = cell_genes.loc[lambda x: x['Paper'] == 'Lake', :]
        elif include == 'lake_plus_glia':
            cell_genes = cell_genes.loc[lambda x: (x['Paper'] == 'Lake') | (~x['Class'].str.contains('Neuro')), :]
        elif include == 'not_lake':
            cell_genes = cell_genes.loc[lambda x: x['Paper'] != 'Lake', :]

        if subtype:
            cell_genes = (cell_genes
            .assign(Class = lambda x: [t if 'Neuro' in c else c for c,t in zip(x['Class'], x['Type'])])
             )
            
        import re
        if combine_layers:
            cell_genes = (cell_genes
            .assign(Class = lambda x: [re.sub('(a|b|c|d|e)', '', c) if bool(re.search('(Ex|In)', c)) else c for c in x['Class']])
             )
            
        if combine_ex_in:
            cell_genes = (cell_genes
             .replace({'Class':{'Neuro-Ex':'Neuro', 'Neuro-In':'Neuro'}})
            )
        else:
            cell_genes = cell_genes.query("Class != 'Neuro'")
            
        cell_genes = (cell_genes
         .melt(id_vars=['Type', 'Paper', 'Cluster', 'Class'], value_name='gene')
         .loc[lambda x: x['gene'].notnull(), ['Class', 'gene']]
         .rename({'Class':'label'}, axis=1)
         .drop_duplicates()
         .loc[lambda x: ~np.isin(x['label'], ['Per'])]
         .sort_values(['label', 'gene'])
         .replace({'gene':replace_dict})
         # .groupby('Class').apply(lambda x: x.sample(frac=.1))
        )

        if add_synapses:
            synapse_genes = get_synapse_genes()
            # cell_genes = cell_genes.loc[lambda x: ~np.isin(x['gene'], synapse_genes['gene'])]
            cell_genes = pd.concat([cell_genes, synapse_genes])
        
    return cell_genes


def get_synapse_genes():
    synapse_genes = pd.read_csv("../data/synaptome_all.csv", usecols=['gene_symbol']).dropna().iloc[:1886]
    synapse_genes = pd.DataFrame({
            'label':'Synapses', 'gene':synapse_genes['gene_symbol']
        }).replace({'gene':replace_dict})
    return synapse_genes


def get_cell_genes_weighted(which=None, normalize=True):
    """
    Read genes table with weights
    """
    
    lake_ex = pd.read_csv("../data/lake_ex.csv")
    lake_in = pd.read_csv("../data/lake_in.csv")
    
    def fisher_norm(y):
        return (np.exp(y) - np.exp(-y))/(np.exp(y) + np.exp(-y))
    
    def clean_lake(df, normalize = True):
        df = (df
         .rename({'cluster':'label', 'Gene':'gene'}, axis=1)
         .melt(id_vars=['label', 'gene'], var_name='which', value_name='weight')
         .loc[lambda x: x['label'] == x['which']].drop('which', axis=1)
         .dropna()
        )
        
        if normalize:
            df = (df
                 .set_index('gene')
                 # .assign(weight = lambda x: x.groupby('label').transform(lambda y: (y-np.mean(y))/np.std(y)))
                 .assign(weight = lambda x: x.groupby('label').transform(lambda y: fisher_norm(y)))
                 .fillna(1)
                 .reset_index()
                 )
        return df
    
    lake_ex = clean_lake(lake_ex, normalize=normalize)
    lake_in = clean_lake(lake_in, normalize=normalize)
    
    if which=='ex':
        return lake_ex
    elif which=='in':
        return lake_in
    else:
        return pd.concat([lake_ex, lake_in])


def get_layer_genes(which='maynard', add_hse_genes=False):
    he_layers = (pd.read_csv("../data/he_layers.csv")
                .loc[:,['Gene symbol', 'Layer marker in human', 'Log2FC to other layers in human']]
                .set_axis(['gene', 'label', 'log2FC'], axis=1)
                .dropna()
                .loc[lambda x: x['log2FC']>1]
                .sort_values('label')
                .drop(['log2FC'],axis=1)
                .replace({'gene':replace_dict})
                )
    
    maynard_data = pd.read_csv("../data/maynard_layers.csv")

    maynard_tstat = (maynard_data
                    .loc[:,['gene', 't_stat_Layer1','t_stat_Layer2','t_stat_Layer3','t_stat_Layer4','t_stat_Layer5','t_stat_Layer6', 't_stat_WM']]
                    .set_index('gene')
                    .set_axis(['L1','L2','L3','L4','L5','L6', 'WM'], axis=1).reset_index()
                    .melt(id_vars='gene', var_name='label', value_name='tstat')
                    .set_index(['gene', 'label'])
                    )

    maynard_layers = (maynard_data
                        .loc[:,['gene', 'fdr_Layer1','fdr_Layer2','fdr_Layer3','fdr_Layer4','fdr_Layer5','fdr_Layer6', 'fdr_WM']]
                        .set_index('gene')
                        .set_axis(['L1','L2','L3','L4','L5','L6', 'WM'], axis=1).reset_index()
                        .melt(id_vars='gene', var_name='label', value_name='fdr')
                        .set_index(['gene', 'label'])
                        .join(maynard_tstat)
                        .loc[lambda x: (x['fdr']<0.05) & (x['tstat']>0)]
                        .reset_index()
                        .sort_values('label')
                        .drop(['tstat', 'fdr'],axis=1)
                        )
    
    he_maynard_layers = pd.concat([
                            he_layers, #.replace({'L6':'L6/WM'}), 
                            maynard_layers.replace({'L6':'L6', 'WM':'L6'})
                        ]).drop_duplicates() 
    

    if which=='he':
        layer_genes = he_layers
    elif which=='maynard':
        layer_genes = maynard_layers
    elif which=='both':
        layer_genes = he_maynard_layers

    hse_genes = ['BEND5',
                'C1QL2',
                'CACNA1E',
                'COL24A1',
                'COL6A1',
                'CRYM',
                'KCNC3',
                'KCNH4',
                'LGALS1',
                'MFGE8',
                'NEFH',
                'PRSS12',
                'SCN3B',
                'SCN4B',
                'SNCG',
                'SV2C',
                'SYT2',
                'TPBG',
                'VAMP1'
                ]
    hse_genes = pd.DataFrame({'gene':hse_genes, 'label':'HSE'})

    if add_hse_genes:
        layer_genes = pd.concat([hse_genes, layer_genes])

    return layer_genes


def get_intelligence_gwas_genes():
    lee2018 = (pd.read_csv("../data/gwas/lee2018_tableS7.csv", header=1, usecols=['Gene symbol'])
               .dropna().rename({'Gene symbol':'gene'}, axis=1).assign(label='lee2018'))
    # jansen2020 = (pd.read_csv("../data/gwas/jansen2020_tableS13.csv", header=2, usecols=['HUGO Symbol'])
    #               .dropna().rename({'HUGO Symbol':'gene'}, axis=1).assign(label='jansen2020'))
    ## UPDATE THIS
    gwas_intelligence = (pd.read_csv("../data/gwas/intelligence.csv", header=0)
        .drop('sniekers2017', axis=1)
        .melt(var_name='label',value_name='gene').dropna()
    )

    rename_dict = {
        # 'lee2018':'Lee 2018',
        # 'sniekers2017':'Sniekers 2017',
        'davies2018':'Davies 2018',
        'savage2018':'Savage 2018',
        'hill2019':'Hill 2019',
        'hatoum2022':'Hatoum 2022',
        # 'jansen2020':'Jansen 2022',
    }

    gwas_intelligence = (
        # pd.concat([gwas_intelligence, lee2018])
        gwas_intelligence
        .replace({'label': rename_dict})
        .assign(label = lambda x: pd.Categorical(x['label'], ordered=True, categories=rename_dict.values()))
    )
    return gwas_intelligence

# def get_structure_gwas_genes():
    # gwas_thickness = (pd.read_csv("../data/gwas/warrier2022_tableS23.csv", index_col=0, header=1).reset_index()
    #              .rename({'Phenotype':'label', 'Gene Symbol':'gene'}, axis=1)
    #              .loc[lambda x: x['label']=='Cortical thickness', ['label','gene']]
    #              )

    # gwas_wm = (pd.read_csv("../data/gwas/sha2023_tableS18.csv", header=1)
    #  .dropna()['Gene name'].rename('gene').to_frame()
    #  .assign(label='WM connectivity')
    # )
    # return pd.concat([gwas_thickness, gwas_wm])
