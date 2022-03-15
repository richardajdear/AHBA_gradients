"""
Version class for analyzing gradients of AHBA
"""

import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSCanonical, PLSRegression
from sklearn.preprocessing import StandardScaler
from brainspace.gradient import GradientMaps
from processing_helpers import *


class gradientVersion():
    
    def __init__(self, n_components=5, approach='pca', sparsity=0, kernel=None, gamma=None, **kwargs):
        """
        Initialize
        """
        self.n_components = n_components
        self.approach = approach
        self.kernel = kernel
        self.sparsity=sparsity
        self.gamma=gamma
        self.kw_emb = kwargs
        self.gradients = GradientMaps(n_components=n_components, approach=approach, kernel=kernel)
        
        self.expression = None
        self.scores = None
        self.var = None

        
    def fit(self, expression, scale=False, message=True, 
            data_dir = "~/rds/rds-cam-psych-transc-Pb9UGUlrwWc/Cam_LIBD/AHBA_data/abagen-data/expression/"):
        """
        Fit to data
        """
        # If expression is given as a string name, read the file
        if isinstance(expression, str):
            X = pd.read_csv(data_dir + expression + '.csv', index_col=0)
        else:
            X = expression
            expression = ''
            
        # Clean data: drop regions with all nulls, and genes with any nulls
        X = X.dropna(axis=0, how='all').dropna(axis=1, how='any')
        
        # Fit gradients
        self.gradients.fit(X.values, sparsity=self.sparsity, gamma=self.gamma, **self.kw_emb)
        
        # Store outputs
        self.expression = X
        self.scores = pd.DataFrame(self.gradients.gradients_, index=X.index)
        self.var = self.gradients.lambdas_
        self.weights = self.fit_weights()
    
        if message:
            print(f"New gradients version: method={self.approach}, kernel={self.kernel}, data={expression}")
        
        return self
    
    def clean_scores(self, flips = [1]):
        """
        Normalize scores, flip as needed, add labels x
        """
        flips = [-1 if i in flips else 1 for i in range(5)]
        scores = self.scores * flips
        
        if self.scores.shape[0]>=120:
            labels = get_labels_hcp()
        elif self.scores.shape[0]<=34:
            labels = get_labels_dk()
        else:
            labels = get_labels_dx()
        
        scores = (scores
                  .apply(lambda x: (x-np.mean(x))/np.std(x))
                  .join(labels)
                 )
        return scores
        
    
    def fit_weights(self, expression=None, independent=True, sort=False, save_name=None, overwrite=True, flips = [0,2], normalize=False):
        """
        Get gene weights from PLS
        Use other expression matrix if provided, otherwise use self
        If independent==True, fit each component separately
        Flip weight gradients as needed
        """
        # Use own expression matrix as X, or use matching regions in another expression matrix
        if expression is None:
            X = self.expression
        else:
            X = expression.dropna(axis=0, how='all').dropna(axis=1, how='any')
            # Make sure X will match onto Y
            X = X.loc[set(X.index).intersection(self.scores.index), :]
            
        # Fit each component independently
        if independent:
            pls_weights = np.zeros((X.shape[1], 5))
            for i in range(5):
                Y = self.scores.loc[X.index, i]
                pls_weights[:,i] = PLSCanonical(n_components=1).fit(X,Y).x_weights_.squeeze()
        # Or fit all components together
        else:
            Y = self.scores.loc[X.index, :]
            pls_weights = PLSCanonical(n_components=5).fit(X,Y).x_weights_
        
        # Normalize
        if normalize:
            pls_weights = StandardScaler().fit_transform(pls_weights)
        
        # Make dataframe and flip as needed
        flips = [-1 if i in flips else 1 for i in range(5)]
        pls_weights = pd.DataFrame(pls_weights, index=X.columns) * flips
        
        # Save to self
        if overwrite:
            self.weights = pls_weights
        
        # Output sorted lists, or unsorted dataframe
        if sort:
            if save_name is not None:
                self.sort_weights(pls_weights).to_csv("../outputs/" + save_name + ".csv")
            return self.sort_weights(pls_weights)
        else:
            return pls_weights
        
        
        
    def sort_weights(self, gene_weights):
        """
        Return genes with each column independently sorted by weight
        """
        gene_ranks = {}    
        for g in range(gene_weights.shape[1]):
            gene_ranks[g] = (
                gene_weights.iloc[:,g]
                .sort_values(ascending=False)
                .reset_index()
                .set_axis(['gene', 'weight'], axis=1)
            )
        
        return pd.concat(gene_ranks, axis=1)

                
    def correlate(self, a,b):
        n = a.shape[1]
        return pd.concat([a,b],axis=1).corr().iloc[:n,n:]
                
        
    def correlate_genes(self, expression, return_ranks=False):
        """
        Get pseudo gene weights by correlating gene expression
        with the scores map
        """
        
        gene_corrs = self.correlate(self.scores, expression)
    
        if return_ranks:
            return self.sort_genes(gene_corrs)
        else:
            return gene_corrs
        
        
    def match_components(self, df_corr):
        """
        Matching logic for a correlation df
        """
        _matches = [None]*5
        _corrs = [None]*5

        for i in range(5):
            # Find columns and rows already matched
            xs = [m[0] for m in _matches if m != None]
            ys = [m[1] for m in _matches if m != None]
            # Take correlation matrix and remove columns and rows already matched
            df_remain = df_corr.drop(xs,axis=0).drop(ys,axis=1)
            # Find the next highest correlation across the remaining values
            new_match = df_remain.stack().index[np.argmax(np.abs(df_remain.values))]
            # Add the new match into the list of matches
            _matches[i] = new_match
            # Add the new correlation into the list of correlations
            _corrs[i] = df_corr.loc[new_match]

        return _matches, _corrs


    def match_and_sort(self, other, df_corr):
        """
        Do correlation matching for genes
        """
        matches, corrs = self.match_components(df_corr)

        xs = [m[0] for m in matches]
        ys = [m[1] for m in matches]
        x_vars = [self.var[i] for i in xs]
        y_vars = [other.var[i] for i in ys]
        mean_vars = [(x+y)/2 for x,y in zip(x_vars, y_vars)]
        sort_idx = np.argsort(mean_vars)[::-1]

        vars_sort = [mean_vars[i] for i in sort_idx]
        corrs_sort = [corrs[i] for i in sort_idx]
        matches_sort = [matches[i] for i in sort_idx]

        return (
            pd.DataFrame(
                [matches_sort, corrs_sort, vars_sort]).T
            .set_axis(['match', 'corr', 'var_expl'], axis=1)
            .astype(dtype={'corr':float, 'var_expl':float})
        )

    def corr_scores(self, other, match=False):
        """
        Correlate region scores with another version, with matching
        """

        df_corr = pd.concat([self.scores, other.scores],axis=1).corr().iloc[:5,5:]
        
        if match:
            return self.match_and_sort(other, df_corr)
        else:
            return df_corr
    
    def corr_weights(self, other, match=False):
        """
        Correlate gene weights with another version
        Optionally with matching
        """
  
        df_corr = pd.concat([self.weights, other.weights],axis=1).corr().iloc[:5,5:]
        
        if match:
            return self.match_and_sort(other, df_corr)
        else:
            return df_corr
    
    

    
    
#     def var_test(self, p=10, k=5):
#         """
#         Do permutation test against variance explained
#         """
#         self.var_perm = np.zeros([p,k])
#         for i in range(p):
#             _df = self.expression.copy().apply(np.random.permutation, axis=0)
#             _var = PCA(n_components=k).fit(_df).explained_variance_ratio_
#             self.var_perm[i,:] = _var
#         print(f'Computed var explained over {p} permutations')
    
    
#     def noise_test(self, reps=10, betas=np.around(np.linspace(.25,5,20),2)):
#         """
#         Do noise test
#         """
#         mu = self.expression.values.mean()
#         sigma = self.expression.values.std()
        
#         # Make arrays for each metric for each PC at each beta
#         var_expl = np.zeros((len(betas), 5))
#         gene_corrs = np.zeros((len(betas), 5))
#         score_corrs = np.zeros((len(betas), 5))
    
        
#         for i, beta in enumerate(betas):
#             _var_expl = np.zeros((reps, 5))
#             _gene_corrs = np.zeros((reps, 5))
#             _score_corrs = np.zeros((reps, 5))
#             for rep in range(reps):
#                 # Add noise and do PCA
#                 expression_noisy = (self.expression + 
#                                     beta * np.random.normal(mu, sigma, self.expression.shape))
#                 pca_rep = PCA(n_components=5).fit(expression_noisy)
#                 coefs_rep = pd.DataFrame(pca_rep.components_, columns=expression_noisy.columns)
#                 scores_rep = pd.DataFrame(pca_rep.transform(expression_noisy), index=expression_noisy.index)
                
#                 # Get metrics for this rep and save into array
#                 _var_expl[rep, :] = pca_rep.explained_variance_ratio_
#                 _gene_corrs[rep, :] = self.coefs.T.corrwith(coefs_rep.T)
#                 _score_corrs[rep, :] = self.scores.corrwith(scores_rep)
                
#             # Take abs mean of all the reps for this beta and save into array
#             var_expl[i, :] = _var_expl.mean(axis=0)
#             gene_corrs[i, :] = abs(_gene_corrs).mean(axis=0)
#             score_corrs[i, :] = abs(_score_corrs).mean(axis=0)

#         self.noise_var_expl = pd.DataFrame(var_expl, index=betas)
#         self.noise_gene_corrs = pd.DataFrame(gene_corrs, index=betas)
#         self.noise_score_corrs = pd.DataFrame(score_corrs, index=betas)
#         print(f'Computed noise metrics for betas = {betas}, abs mean over {reps} reps')

 

#     def missing_test(self, reps=40, missing=[.1,.2,.3,.4,.5], match=True):
#         """
#         Do missing regions test
#         """
        
#         # Make arrays for each metric for each PC at each beta
#         var_expl = np.zeros((len(missing), 5))
#         gene_corrs = np.zeros((len(missing), 5))
#         score_corrs = np.zeros((len(missing), 5))
        
#         for i, miss in enumerate(missing):
#             _var_expl = np.zeros((reps, 5))
#             _gene_corrs = np.zeros((reps, 5))
#             _score_corrs = np.zeros((reps, 5))
#             for rep in range(reps):
#                 # Add noise and do PCA
#                 expression_missing = (self.expression.sample(frac=1-miss, replace=False))
#                 pca_rep = pcaVersion(expression_missing, message=False)
                
#                 # Get metrics for this rep and save into array
#                 _var_expl[rep, :] = pca_rep.var
#                 if match:
#                     _gene_corrs[rep, :] = self.corr_coefs(pca_rep, match=True)['corr'].values
#                     _score_corrs[rep, :] = self.corr_scores(pca_rep, match=True)['corr'].values
#                 else:
#                     _gene_corrs[rep, :] = self.corr_coefs(pca_rep).pipe(np.diag)
#                     _score_corrs[rep, :] = self.corr_scores(pca_rep).pipe(np.diag)
                
#             # Take abs mean of all the reps for this beta and save into array
#             var_expl[i, :] = _var_expl.mean(axis=0)
#             gene_corrs[i, :] = abs(_gene_corrs).mean(axis=0)
#             score_corrs[i, :] = abs(_score_corrs).mean(axis=0)
        
#         if match:
#             self.missing_var_expl = pd.DataFrame(var_expl, index=missing)
#             self.missing_gene_corrs = pd.DataFrame(gene_corrs, index=missing)
#             self.missing_score_corrs = pd.DataFrame(score_corrs, index=missing)
#         elif ~match:
#             self.missing_nomatch_var_expl = pd.DataFrame(var_expl, index=missing)
#             self.missing_nomatch_gene_corrs = pd.DataFrame(gene_corrs, index=missing)
#             self.missing_nomatch_score_corrs = pd.DataFrame(score_corrs, index=missing)
            
#         print(f'Computed missing metrics after dropping {missing} regions, match={match}, abs mean over {reps} reps')

    
#     def bootstrap_test(self, q=100):
#         """
#         Do bootstrap test
#         """
        
#         X = self.expression
#         _coef_corrs = {'match': np.zeros((q, 5)), 
#                        'nomatch': np.zeros((q, 5))}
#         _score_corrs = {'match': np.zeros((q, 5)), 
#                        'nomatch': np.zeros((q, 5))}
#         for i in range(q):
#             _exp = X.sample(frac=1,replace=True)
#             _pca = pcaVersion(_exp, message=False)
#             _pca.scores = _pca.scores.groupby('label').mean() # aggregate over duplicate regions
#             _coef_corrs['match'][i,:] = _pca.corr_coefs(self, match=True)['corr'].values
#             _score_corrs['match'][i,:] = _pca.corr_scores(self, match=True)['corr'].values
#             _coef_corrs['nomatch'][i,:] = _pca.corr_coefs(self, match=False).pipe(np.diag)
#             _score_corrs['nomatch'][i,:] = _pca.corr_scores(self, match=False).pipe(np.diag)

#         self.boot_test_coefs = {k:pd.DataFrame(v) for k,v in _coef_corrs.items()}
#         self.boot_test_scores = {k:pd.DataFrame(v) for k,v in _score_corrs.items()}
#         print('Bootstrap test done')
    
    
    
#     def bootstrap_coefs(self, q=1000):
#         """
#         Get bootstrapped coefficients
#         """
#         X = self.expression
#         _coefs = np.zeros((5, X.shape[1], q))
#         for i in range(q):
#             _exp = X.sample(frac=1,replace=True)
#             _pca = PCA(n_components=5).fit(_exp)
#             _coefs[:,:,i] = _pca.components_
# #             if i%100==0: 
# #                 print(f'Bootstrapped {i} times')

#         coefs_boot = np.mean(_coefs, axis=2)
#         coefs_boot_norm = np.apply_along_axis(lambda x: np.mean(x/np.std(x)), axis=2, arr=_coefs)

#         self.coefs_boot = pd.DataFrame(coefs_boot, columns=X.columns)
#         self.coefs_boot_norm = pd.DataFrame(coefs_boot_norm, columns=X.columns)
#         print('Bootstrap done!')
    
    
    
### LEGACY
    
    
#     # Project two sets of coefs into roi-space of this pcaVersion and get correlation
#     # Coefs could be from this or other pcaVersions
#     # Rotates scores2 onto scores1
    def score_and_corr(self, version1, version2, procrustes=False):
        
        scores1 = self.score_from(version1)
        
        if not procrustes:
            scores2 = self.score_from(version2)
        elif procrustes:
            scores2 = self.expression @ (version1.rotate_other_coefs(version2).T)
        
        return (
            pd.concat([scores1, scores2],axis=1)
            .corr().iloc[:5,5:]
        )
    
    

    def varimax(self):
        """
        Apply varimax to self and return rotation matrix R
        """
        L,R = rotate_factors(self.loadings.values, method='varimax')
        p,k = L.shape
        return R
    
    def quartimax(self, how='loadings'):
        """
        Apply quartimax to self and return rotation matrix R
        """
        if how=='loadings':
            L,R = rotate_factors(self.loadings.values, method='quartimax')
        elif how=='scores':
            L,R = rotate_factors(self.U.values, method='quartimax')
        p,k = L.shape
        return R

    
    def procrustes(self, other, varimax=False, gamma=1, q=100):
        """
        Apply procrustes rotation to self to match target other and return rotation matrix
        """
        if not varimax:
            L1 = self.loadings
            L2 = other.loadings
        elif varimax:
            L1 = self.loadings @ self.varimax(gamma=gamma, q=q)
            L2 = other.loadings @ other.varimax(gamma=gamma, q=q)
        
        p,k = L1.shape
        R, _ = orthogonal_procrustes(L1, L2)
#         coefs_rotated = (self.coefs.T @ R_dim5).set_axis(list(range(1,6)),axis=1)
        return R

    
    def corr_loadings(self, other, varimax=False, quartimax=False, gamma=1, q=100, procrustes=False, target=None):
        """
        Get correlation df of top 5 components with another version x
        Optionally apply varimax first
        Optionally Procrustes-rotate self onto other or onto target, with or without varimax
        """
        
        if varimax:
            L1 = self.loadings @ self.varimax(gamma=gamma, q=q)
            L2 = other.loadings @ other.varimax(gamma=gamma, q=q)
        elif quartimax:
            L1 = self.loadings @ self.quartimax()
            L2 = other.loadings @ other.quartimax()
        else:
            L1 = self.loadings
            L2 = other.loadings
        
        if procrustes and target is None:
            L1 = L1 @ self.procrustes(other, varimax=varimax, gamma=gamma, q=q)
        elif procrustes and target is not None:
            L1 = L1 @ self.procrustes(target, varimax=varimax, gamma=gamma, q=q)
            L2 = L2 @ other.procrustes(target, varimax=varimax, gamma=gamma, q=q)
        
        return (
            pd.concat([L1.iloc[:,:5], L2.iloc[:,:5]],axis=1) # only look at top 5 components
            .corr().iloc[:5,5:]
        )
    
    

    
    


