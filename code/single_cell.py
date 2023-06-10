import numpy as np, pandas as pd

# Single cell data analysis

def get_single_cell_matching_genes(weights, single_cell_path="../data/allen_single_cell_matrix.csv"):
    genes_to_match = set(weights.index)
    all_columns = pd.read_csv(single_cell_path, nrows=1).columns
    matching_column_indices = [i for i, e in enumerate(all_columns) if e in genes_to_match]
    weights_matched = weights.loc[list(set(all_columns).intersection(genes_to_match)),:]
    return weights_matched, matching_column_indices


def get_weights_posneg(weights):
    weights_posneg = pd.concat([
        weights.applymap(lambda x: np.where(x<0, -x, 0)).set_axis([f'C{i+1}-' for i in range(3)], axis=1),
        weights.applymap(lambda x: np.where(x>0, x, 0)).set_axis([f'C{i+1}+' for i in range(3)], axis=1)
    ], axis=1)
    return weights_posneg

def project_single_cell_posneg(weights_matched_posneg, matching_column_indices, 
                               single_cell_path = "../data/allen_single_cell_matrix.csv",
                               save_path = "../outputs/allen_single_cell_C123_posneg.csv"):
    sc_projected_posneg = (
        pd.read_csv(single_cell_path, header=0, usecols=matching_column_indices)
        .add(1)
        .pipe(np.log10)
        .pipe(lambda x: x/x.mean(axis=1).values[:,None])
        .dot(weights_matched_posneg)
    )
    sc_projected_posneg.to_csv(save_path)
    return sc_projected_posneg
    

def format_single_cell_projected_for_plot(sc_projected_posneg = "../outputs/allen_single_cell_C123_posneg.csv",
                                          metadata_path = "../data/allen_single_cell_metadata.csv"):
    """
    Format Allen Single Cell data for positive vs negative plot
    NB: must first project gene data to C123
    """
    if isinstance(sc_projected_posneg, str):
        sc_projected_posneg = pd.read_csv(sc_projected_posneg, index_col=0)

    metadata = (pd.read_csv(metadata_path)
                .assign(cell_type = lambda x: np.where(x['class_label']=='Non-neuronal', x['subclass_label'], 
                                 x['class_label'].map({'Glutamatergic':'Neuron-Ex','GABAergic':'Neuron-In'})))
                .assign(cell_type = lambda x: pd.Categorical(x['cell_type'], ordered=True,
                                           categories=['Neuron-Ex','Neuron-In','Astrocyte','Endothelial','Microglia','Oligodendrocyte','OPC']))
                .assign(layer = lambda x: x['cortical_layer_label'].str[:2])
    )
                                
    sc_projected_pos = sc_projected_posneg.iloc[:,:3].set_axis(['C1','C2','C3'],axis=1).stack().rename('negative')
    sc_projected_neg = sc_projected_posneg.iloc[:,3:].set_axis(['C1','C2','C3'],axis=1).stack().rename('positive')
    sc_projected_melt = sc_projected_pos.to_frame().join(sc_projected_neg).reset_index(1).rename({'level_1':'C'}, axis=1)

    cell_rename = {
        'Neuron-Ex':'N-Ex',
        'Neuron-In':'N-In',
        'Astrocyte':'Astro',
        'Microglia':'Micro',
        'Endothelial':'Endo',
        'Oligodendrocyte':'Oligo',
        'OPC':'OPC'
    }

    sc_projected_posneg_plot = (
        metadata.loc[:,['class_label','subclass_label','region_label','layer','cell_type','outlier_call']]
        .join(sc_projected_melt)
        .query("~outlier_call")
        # .query("cell_type not in ['Pericyte', 'VLMC']")
        .sort_values("cell_type")
        .dropna()
        .replace({'cell_type':cell_rename})
        # .reset_index(drop=True)
    )
    return sc_projected_posneg_plot
