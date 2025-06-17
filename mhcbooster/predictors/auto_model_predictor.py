import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from adjustText import adjust_text
from sklearn.decomposition import TruncatedSVD


def _pred_best_combination(score_matrix, report_directory, predictor_type=None):

    predictors = np.array([column.replace('_log_rt_error', '').replace('_entropy_score','') for column in score_matrix.columns])
    best_predictors = []

    svd = TruncatedSVD(n_components=score_matrix.shape[1])
    svd.fit(score_matrix)

    feature_importance_matrix = svd.singular_values_.reshape(-1,1) * svd.components_
    feature_importance_matrix[0] *= -1 if np.max(feature_importance_matrix[0]) < 0 else 1

    ### Noise filtering
    dim_importance_ratio = (svd.singular_values_ ** 2) / np.sum(svd.singular_values_ ** 2)
    cumulative_importance_variance = np.cumsum(dim_importance_ratio)
    n_components = np.argmax(cumulative_importance_variance >= 0.99) + 1
    best_model_index = np.argmax(feature_importance_matrix[0])
    best_predictors.append(predictors[best_model_index])

    ### Model filtering
    if n_components > 1:
        best_score_dim0 = feature_importance_matrix[0][best_model_index]
        avail_indices = np.where(feature_importance_matrix[0] >= 0.95 * best_score_dim0)[0]

        ### Find best models
        if len(avail_indices) > 1:
            truncated_importance_matrix = feature_importance_matrix[1:n_components, avail_indices]
            truncated_importance_matrix = truncated_importance_matrix - feature_importance_matrix[1:n_components,best_model_index].reshape(-1,1)
            distance =np.linalg.norm(truncated_importance_matrix, axis=0)
            furthest_index = np.argmax(distance)
            best_predictors.append(predictors[avail_indices][furthest_index])

    ### Draw distribution
    colors = ['red' if p in best_predictors else 'blue' for p in predictors]
    plt.scatter(feature_importance_matrix[0], feature_importance_matrix[1], c=colors)
    texts = []
    for xi, yi, label in zip(feature_importance_matrix[0], feature_importance_matrix[1], predictors):
        texts.append(plt.text(xi, yi, label, fontsize=9, ha='left', va='top'))
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray'))
    plt.savefig(report_directory / f'{predictor_type}_model_representativeness.png')
    plt.close()
    # plt.show()

    print(f"Best {predictor_type} predictors: ", best_predictors)
    return best_predictors

def predict_best_combination(feature_matrix, report_directory):
    rt_scores = feature_matrix[[col for col in feature_matrix.columns if 'log_rt_error' in col and 'Chronologer' not in col]]
    ms2_scores = feature_matrix[[col for col in feature_matrix.columns if 'entropy_score' in col]]

    rt_predictors = _pred_best_combination(rt_scores, report_directory, predictor_type='RT')
    ms2_predictors = _pred_best_combination(ms2_scores, report_directory, predictor_type='MS2')

    return rt_predictors, ms2_predictors

if __name__ == '__main__':
    feature_matrix = pd.read_csv('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/best/JY_Class1_25M_DDA_60min_Slot1-12_1_552_MHCBooster/all_features.tsv', sep='\t')
    feature_matrix = pd.read_csv('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M_rerun/msfragger/mhcbooster_newauto/JY_Class1_25M_DDA_60min_Slot1-12_1_552_MHCBooster/all_features.tsv', sep='\t')
    predict_best_combination(feature_matrix, None)
