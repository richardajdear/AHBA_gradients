"""
Version class for analyzing PCA of AHBA donor triplets
Builds off pcaVersion class
"""

import numpy as np, pandas as pd
from itertools import combinations

import abagen
from pcaVersion import pcaVersion


class tripletsVersion():
    
    def __init__(self, expression_donors, k=5):
        donors = ['9861', '10021', '12876', '14380', '15496', '15697']
        triplets_list = [list(x) for x in list(combinations(range(6), 3))]
        triplets_names = [''.join(map(str,x)) for x in triplets_list]
        self.triplets_dict = dict(zip(triplets_names, triplets_list))
        self.disjoint_triplets = [list(x) for x in combinations(triplets_names,2) 
                             if not any(set(list(x[0])).intersection(set(list(x[1]))))]
        
        self.expression_donors = expression_donors
        self.triplets_pca = self.run_pca(k=k)
        
        print("New PCA triplets")
        
    
    def run_pca(self, k=5):
        """
        Run PCA for each triplet x
        """
        triplets_pca = {}
        for name, triplet in self.triplets_dict.items():
            # Get the data from the three donors
            expression_triplet = [self.expression_donors[i] for i in triplet]
            # Combine the data
            data = pd.concat(expression_triplet).groupby(level=0).mean()
            # Do the PCA
            #print("Running PCA for donors %s" % name)
            triplets_pca[name] = pcaVersion(data, k=k, message=False)
            
        return triplets_pca
            
            
    def disjoint_corrs(self, how='loadings', procrustes_target=None, quartimax=False, match_components=False, DS=False):
        self.corrs = {}
        
        for pair in self.disjoint_triplets:
            name = '-'.join(pair)
            if not DS:
                pca1 = self.triplets_pca[pair[0]]
                pca2 = self.triplets_pca[pair[1]]
            elif DS:
                pca1 = self.triplets_pca_DS[pair[0]]
                pca2 = self.triplets_pca_DS[pair[1]]
            
            if how=='loadings':
                if procrustes_target is not None:
                    df_corr = pca1.corr_loadings(pca2, procrustes=True, target=procrustes_target, quartimax=quartimax)
                else:
                    df_corr = pca1.corr_loadings(pca2, quartimax=quartimax)
            elif how=='scores':
                if procrustes_target is not None:
                    df_corr = pca1.corr_scores(pca2, procrustes=True, target=procrustes_target)
                else:
                    df_corr = pca1.corr_scores(pca2)
                
            if match_components:
                self.corrs[name] = pca1.corr_match_sort(pca2, df_corr)['corr']
            else:
                self.corrs[name] = df_corr.pipe(np.diag)
            
        return pd.DataFrame(self.corrs, index=list(range(1,6)))
    
    
    def filter_differential_stability(self, threshold=.1, k=5):
        self.triplets_pca_DS = {}
        for name, triplet in self.triplets_dict.items():
            # Get the data from the three donors
            expression_triplet = [self.expression_donors[i] for i in triplet]
            # Filter genes by stability in this triplet
            expression_triplet_DS = abagen.correct.keep_stable_genes(expression_triplet, threshold=threshold)
            # Combine the data
            data = pd.concat(expression_triplet_DS).groupby(level=0).mean()
            # Do the PCA
            #print("Running PCA for donors %s" % name)
            self.triplets_pca_DS[name] = pcaVersion(data, k=k, message=False)
        print(f"Computed differential stability PCAs, threshold={threshold}")