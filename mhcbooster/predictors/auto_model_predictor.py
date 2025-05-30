import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD


def _pred_best_combination(score_matrix, predictor_type=None, max_predictors=4):

    predictors = np.array([column.replace('_log_rt_error', '').replace('_entropy_score','') for column in score_matrix.columns])

    ### First round decomposition: noise filtering
    svd = TruncatedSVD(n_components=score_matrix.shape[1])
    svd.fit(score_matrix)

    # Calculate contribution of each principle component
    explained_variance_ratio = (svd.singular_values_ ** 2) / np.sum(svd.singular_values_ ** 2)
    cumulative_explained_variance = np.cumsum(explained_variance_ratio)

    # Determine n_components by the smallest number with a cumulative contribution rate of 90%
    n_components = np.argmax(cumulative_explained_variance >= 0.9) + 1
    n_components = min(max(2, n_components), max_predictors)

    ### Second round decomposition: choose the most important scores
    # svd = TruncatedSVD(n_components=n_components)
    # svd.fit(score_matrix)

    # Compute the L2 norm of each column (Euclidean norm) for "importance"
    importance_matrix = np.dot(np.diag(svd.singular_values_[:n_components]), svd.components_[:n_components,])
    importance = np.linalg.norm(importance_matrix, axis=0)
    sort_indices = np.argsort(-importance)
    best_predictors = list(predictors[sort_indices][: n_components])
    print(f"Best {predictor_type} predictors: ", best_predictors)
    return best_predictors


def predict_best_combination(feature_matrix):
    rt_scores = feature_matrix[[col for col in feature_matrix.columns if 'log_rt_error' in col and 'Chronologer' not in col]]
    ms2_scores = feature_matrix[[col for col in feature_matrix.columns if 'entropy_score' in col]]

    rt_predictors = _pred_best_combination(rt_scores, predictor_type='RT')
    ms2_predictors = _pred_best_combination(ms2_scores, predictor_type='MS2')

    return rt_predictors, ms2_predictors

if __name__ == '__main__':
    feature_matrix = pd.read_csv('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/best/JY_Class1_25M_DDA_60min_Slot1-12_1_552_MHCBooster/all_features.tsv', sep='\t')
    predict_best_combination(feature_matrix)