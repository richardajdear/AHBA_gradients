"""
Version class for analyzing gradients of AHBA
"""
from unittest.mock import patch
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSCanonical, PLSRegression
from sklearn.preprocessing import StandardScaler
import brainspace
from brainsmash.mapgen.base import Base
from neuromaps.images import annot_to_gifti
from neuromaps.nulls.spins import parcels_to_vertices, vertices_to_parcels

from processing_helpers import *


### Monkey Patch compute_affinity function to not zero-out negative values
def compute_affinity_new(x, kernel=None, sparsity=.9, pre_sparsify=True,
                    non_negative=False, gamma=None):

    # run original function with different output
    return brainspace.gradient.kernels.compute_affinity(x, kernel=kernel, 
                sparsity=sparsity, pre_sparsify=pre_sparsify, non_negative=False, gamma=gamma)

# NB: must patch in 'gradient.kernels' namespace, not 'gradient.gradient'
brainspace.gradient.gradient.compute_affinity = compute_affinity_new


class gradientVersion():

    
    def __init__(self, approach='dm', n_components=5, sparsity=0, kernel=None,
                #  marker_genes=['NEFL', 'LGALS1', 'SYT6'], 
                 marker_genes=['NEFL', 'LGALS1', 'RFTN1'], 
                **kwargs):
        """
        Initialize
        """
        self.n_components = n_components
        self.marker_genes = marker_genes
        
        self.approach = approach
        self.sparsity = sparsity

        if approach == 'dm':
            kwargs['alpha'] = kwargs.get('alpha', 1) # set alpha=1 as default, but only for approach = 'dm'
            kernel = 'normalized_angle' if kernel is None else kernel
        self.params = kwargs # set embedding-specific parameters
        self.kernel = kernel

        self.expression = None
        self.scores = None
        self.gradients = brainspace.gradient.gradient.GradientMaps(n_components=n_components, approach=approach, kernel=kernel)    

        
    def fit(self, expression, scale=False, message=True, 
            data_dir = "../data/abagen-data/expression/"):
        """
        Fit to data
        Marker genes is a list of genes to define the gradient direction, if gradient n is inversely aligned to marker n, the gradient will be flipped
        """
        # If expression is given as a string name, read the file
        if isinstance(expression, str):
            X = pd.read_csv(data_dir + expression + '.csv', index_col=0)
        else:
            X = expression
            # expression = name # for printing output

        # Clean data: drop regions with all nulls, and genes with any nulls
        X = X.dropna(axis=0, how='all').dropna(axis=1, how='any')
        # Optionally z-score expression
        if scale:
            X = X.apply(lambda x: (x-np.mean(x))/np.std(x))
        self.expression = X

        # Fit gradients
        self.gradients.fit(X.values, sparsity=self.sparsity, **self.params)
        
        # Align gradients to marker genes
        scores = pd.DataFrame(self.gradients.gradients_, index=X.index)
        for i, marker in enumerate(self.marker_genes):
            r_ = X.loc[:,marker].corr(scores.loc[:,i])
            if r_ < 0:
                scores.loc[:,i] *= -1
        
        self.scores = scores
        self.affinity = compute_affinity_new(x=X.values, kernel=self.kernel, sparsity=self.sparsity)
        self.eigenvalues = self.gradients.lambdas_
        self.weights = self.fit_weights()
    
        if message:
            data_name = expression if isinstance(expression, str) else '(data given)'
            print(f"New gradients version: method={self.approach}, kernel={self.kernel}, sparsity={self.sparsity}, data={data_name}")
        
        return self

    
    def clean_scores(self, scores=None, norm=True, n_components=3):
        """
        Normalize G1-3 scores, add labels
        """
        if scores is None:
            scores = self.scores
        
        if scores.shape[0]>=120:
            labels = get_labels_hcp()
        elif scores.shape[0]<=34:
            labels = get_labels_dk()
        else:
            labels = get_labels_dx()
        
        scores = (scores
                  .iloc[:,:n_components]
                  .set_axis(['G'+str(i+1) for i in range(n_components)],axis=1)
                  .rename_axis('id')
                 )
        
        if norm:
            scores = scores.apply(lambda x: (x-np.mean(x))/np.std(x))

        return scores.join(labels)
    
    
    def score_from(self, other, clean=True):
        """
        Score from other weights
        """
        other_weights = other.weights.set_axis(range(other.weights.shape[1]), axis=1)
        gene_intersection = np.intersect1d(self.expression.columns, other_weights.index)
        _expression = self.expression.loc[:, gene_intersection]
        _weights = other_weights.loc[gene_intersection, :]
        
        scores = _expression @ _weights
        
        if clean:
            scores = self.clean_scores(scores=scores)
            
        return scores


    def fill_missing_scores(self, other, clean=True):
        """
        Fill in NA regions in this set of scores using scores from other version
        """
        filled_scores = (self.scores
                         .reindex(range(1,181))
                         .fillna(other.scores)
                         .dropna()
                         )
        
        if clean:
            filled_scores = self.clean_scores(scores=filled_scores)
            
        return filled_scores
    
    def score_in_dk(self, clean=True,
                    hcp_img_path = "../data/parcellations/lh.HCPMMP1.annot",
                    dk_img_path = "../data/parcellations/lh.aparc.annot"
                   ):
        """
        Project HCP scores to a different parcellation using annot files
        """
        hcp_img = annot_to_gifti(hcp_img_path)
        dk_img = annot_to_gifti(dk_img_path)
        
        scores_dk = np.zeros((34,3))
        for i in range(3):
            # Re-index gradient null values with NA
            g_hcp = self.scores[i].reindex(range(1,181)).values
            # Use HCP parcellation image to project HCP data to fsaverage
            g_fsaverage = parcels_to_vertices(g_hcp, hcp_img)
            # Use DK parcellation image to project fsaverage data into DK
            g_dk = vertices_to_parcels(g_fsaverage, dk_img)
            # Add to outputs
            scores_dk[:,i] = g_dk
        
        # Convert to dataframe
        scores_dk = pd.DataFrame.from_records(scores_dk, index=list(range(1,35)))
        
        if clean:
            scores = self.clean_scores(scores=scores_dk)
        
        return scores

    
    def fit_weights(self, sort=False, n_components=3):
        """
        Get gene weights by correlating expression with scores
        """
        x = self.expression.values
        y = self.scores.values
        xv = x - x.mean(axis=0)
        yv = y - y.mean(axis=0)
        xvss = (xv * xv).sum(axis=0)
        yvss = (yv * yv).sum(axis=0)
        result = np.matmul(xv.transpose(), yv) / np.sqrt(np.outer(xvss, yvss))
        # bound the values to -1 to 1 in the event of precision issues
        result = np.maximum(np.minimum(result, 1.0), -1.0)

        weights = pd.DataFrame(result, index=self.expression.columns, 
                               columns=['G'+str(i+1) for i in range(5)]).iloc[:,:n_components]

        # Output sorted lists, or unsorted dataframe
        if sort:
            return self.sort_weights(weights)
        else:
            return weights



    def fit_weights_PLS(self, other_expression=None, independent=True, normalize=False, sort=False, save_name=None, overwrite=True):
        """
        Get gene weights from PLS
        Use other expression matrix if provided, otherwise use self
        If independent==True, fit each component separately
        """
        # Use own expression matrix as X, or use matching regions in another expression matrix
        if other_expression is None:
            X = self.expression
        else:
            X = other_expression.dropna(axis=0, how='all').dropna(axis=1, how='any')
            # Make sure X will match onto regions in Y
            X = X.loc[np.intersect1d(X.index, self.scores.index), :]
        
        # Fit each component independently
        n_components = self.scores.shape[1]
        if independent:
            pls_weights = np.zeros((X.shape[1], n_components))
            for i in range(n_components):
                Y = self.scores.loc[X.index, i]
                pls_weights[:,i] = PLSCanonical(n_components=1).fit(X,Y).x_weights_.squeeze()
                # pls_weights[:,i] = PLSCanonical(n_components=1).fit(X,Y).x_loadings_.squeeze()
        # Or fit all components together
        else:
            Y = self.scores.loc[X.index, :]
            # pls_weights = PLSCanonical(n_components=5).fit(X,Y).x_weights_
            pls_weights = PLSCanonical(n_components=5).fit(X,Y).x_loadings_
        
        # Normalize
        if normalize:
            pls_weights = StandardScaler().fit_transform(pls_weights)
        
        # Make dataframe and align to marker genes
        pls_weights = pd.DataFrame(pls_weights, index=X.columns)
        for i, marker in enumerate(self.marker_genes):
            # If the marker for a component is negative, flip that component
            if pls_weights.loc[marker, i] < 0:
                pls_weights.loc[:,i] *= -1
        
        # Save to self
        if overwrite:
            self.weights = pls_weights
        
        # Output sorted lists, or unsorted dataframe
        if sort:
            if save_name is not None:
                self.sort_weights(pls_weights).to_csv("../outputs/" + save_name + ".csv")
            return self.sort_weights(pls_weights)
        else:
            return pls_weights.set_axis(['G'+str(i+1) for i in range(n_components)], axis=1)


        
        
    def make_null_scores(self, n=10, atlas='hcp', dist_mat=None, save_name=None):
        """
        Generate null maps using brainsmash
        """
        # Choose distance matric that matches parcellation
        if atlas == 'hcp' and dist_mat is None:
            dist_mat=np.loadtxt("../data/LeftParcelGeodesicDistmat.txt")
        elif atlas == 'dk':
            dist_mat=np.loadtxt("../data/LeftParcelGeodesicDistmat_DK.txt")
        # Filter distance matrix to non-null regions
        scores = self.clean_scores()
        inds = [i-1 for i in scores.index]
        dist_mat = dist_mat[inds,:][:,inds]

        null_scores = np.zeros([scores.shape[0], 3, n])

        for m in range(3):
            base_map = scores.iloc[:,m].values
            base = Base(x=base_map, D=dist_mat)
            nulls = base(n)
            null_scores[:,m,:] = nulls.swapaxes(0,1)

        if save_name is not None:
            outfile = "../outputs/permutations/" + save_name + '.npy'
            np.save(outfile, null_scores)
        
        return null_scores
        

    
    def make_null_weights(self, null_scores, save_name = None):
        """
        Generate null weights from PLS on null maps
        """
        n_genes = self.weights.shape[0]
        n = null_scores.shape[2]
        null_weights = np.zeros((n_genes, 3, n))
        
        X = self.expression
        
        for i in range(3):
            # Get index of marker gene to align PLS signs
            marker_ix = X.columns.tolist().index(self.marker_genes[i])
            for m in range(n):
                Y = null_scores[:,i,m]
                nan_mask = np.isnan(Y)
                Y_nan = Y[~nan_mask]
                X_nan = X.values[~nan_mask, :]
                
                _null_weights = PLSCanonical(n_components=1).fit(X_nan,Y_nan).x_weights_.squeeze()
                # Flip output if needed
                if _null_weights[marker_ix] < 0:
                    _null_weights *= -1
                null_weights[:,i,m] = _null_weights
        
        if save_name is not None:
            outfile = "../outputs/permutations/" + save_name + '.npy'
            np.save(outfile, null_weights)
        
        return null_weights
                
        
        
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



    
                
    def correlate(self, a, b):
        # n = a.shape[1]
        # return pd.concat([a,b],axis=1).corr().iloc[:n,n:]

        n = len(a)
        a, b = a.values, b.values
        sums = np.multiply.outer(b.sum(), a.sum())
        stds = np.multiply.outer(b.std(), a.std())
        corr = pd.DataFrame((b.T.dot(a) - sums / n) / stds / n, 
                            b.columns, a.columns)
        return corr

        
                
        
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
        
        
    def match_components(self, df_corr, n_components=5):
        """
        Matching logic for a correlation df
        """
        _matches = [None]*n_components
        _corrs = [None]*n_components

        for i in range(n_components):
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
        x_vars = [self.eigenvalues[i] for i in xs]
        y_vars = [other.eigenvalues[i] for i in ys]
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
#     def score_and_corr(self, version1, version2, procrustes=False):
        
#         scores1 = self.score_from(version1)
        
#         if not procrustes:
#             scores2 = self.score_from(version2)
#         elif procrustes:
#             scores2 = self.expression @ (version1.rotate_other_coefs(version2).T)
        
#         return (
#             pd.concat([scores1, scores2],axis=1)
#             .corr().iloc[:5,5:]
#         )
    
    

#     def varimax(self):
#         """
#         Apply varimax to self and return rotation matrix R
#         """
#         L,R = rotate_factors(self.loadings.values, method='varimax')
#         p,k = L.shape
#         return R
    
#     def quartimax(self, how='loadings'):
#         """
#         Apply quartimax to self and return rotation matrix R
#         """
#         if how=='loadings':
#             L,R = rotate_factors(self.loadings.values, method='quartimax')
#         elif how=='scores':
#             L,R = rotate_factors(self.U.values, method='quartimax')
#         p,k = L.shape
#         return R

    
#     def procrustes(self, other, varimax=False, gamma=1, q=100):
#         """
#         Apply procrustes rotation to self to match target other and return rotation matrix
#         """
#         if not varimax:
#             L1 = self.loadings
#             L2 = other.loadings
#         elif varimax:
#             L1 = self.loadings @ self.varimax(gamma=gamma, q=q)
#             L2 = other.loadings @ other.varimax(gamma=gamma, q=q)
        
#         p,k = L1.shape
#         R, _ = orthogonal_procrustes(L1, L2)
# #         coefs_rotated = (self.coefs.T @ R_dim5).set_axis(list(range(1,6)),axis=1)
#         return R

    
#     def corr_loadings(self, other, varimax=False, quartimax=False, gamma=1, q=100, procrustes=False, target=None):
#         """
#         Get correlation df of top 5 components with another version x
#         Optionally apply varimax first
#         Optionally Procrustes-rotate self onto other or onto target, with or without varimax
#         """
        
#         if varimax:
#             L1 = self.loadings @ self.varimax(gamma=gamma, q=q)
#             L2 = other.loadings @ other.varimax(gamma=gamma, q=q)
#         elif quartimax:
#             L1 = self.loadings @ self.quartimax()
#             L2 = other.loadings @ other.quartimax()
#         else:
#             L1 = self.loadings
#             L2 = other.loadings
        
#         if procrustes and target is None:
#             L1 = L1 @ self.procrustes(other, varimax=varimax, gamma=gamma, q=q)
#         elif procrustes and target is not None:
#             L1 = L1 @ self.procrustes(target, varimax=varimax, gamma=gamma, q=q)
#             L2 = L2 @ other.procrustes(target, varimax=varimax, gamma=gamma, q=q)
        
#         return (
#             pd.concat([L1.iloc[:,:5], L2.iloc[:,:5]],axis=1) # only look at top 5 components
#             .corr().iloc[:5,5:]
#         )
    
    

    
    


