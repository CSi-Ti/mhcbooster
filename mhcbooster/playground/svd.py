import numpy as np
import pandas as pd
from sklearn.decomposition import TruncatedSVD
import matplotlib.pyplot as plt
from mhcbooster.utils.fdr import calculate_qs
from pathlib import Path
from mhcbooster.predictors.auto_model_predictor import predict_best_combination

feature_paths = list(Path('/mnt/d/workspace/mhc-validator-2/experiment/MSV000091456/A375_lowInput_IP/HLA-II/test').rglob('*.features.tsv'))

for i in range(len(feature_paths)):
    feature_matrix = pd.read_csv(feature_paths[0], sep='\t')
    predictions = predict_best_combination(feature_matrix)
    # rt_scores = feature_matrix[[col for col in feature_matrix.columns if 'log_rt_error' in col and 'Chronologer' not in col]]
    # # rt_scores = feature_matrix[[col for col in feature_matrix.columns if 'entropy' in col]]
    # # rt_scores = feature_matrix[[col for col in feature_matrix.columns if 'im_error' in col]]
    #
    # # for col in rt_scores.columns:
    # #     qs = calculate_qs(rt_scores[col], feature_matrix['Label'], higher_better=True)
    # #     print(col, np.sum(qs < 0.01))
    # # if 'lnExpect' in feature_matrix.columns:
    # #     qs = calculate_qs(feature_matrix['lnExpect'].astype(float), feature_matrix['Label'], higher_better=False)
    # # elif 'log10_evalue' in feature_matrix.columns:
    # #     qs = calculate_qs(feature_matrix['log10_evalue'].astype(float), feature_matrix['Label'], higher_better=False)
    # # rt_scores = rt_scores[qs < 0.01]
    #
    # # 执行SVD
    # svd = TruncatedSVD(n_components=rt_scores.shape[1])  # 设置最大分量数为特征数
    # svd.fit(rt_scores)
    #
    # # 计算每个主成分的方差贡献率
    # explained_variance_ratio = (svd.singular_values_ ** 2) / np.sum(svd.singular_values_ ** 2)
    # cumulative_explained_variance = np.cumsum(explained_variance_ratio)
    # print(f"累计方差贡献率: {cumulative_explained_variance}")
    #
    # # 选择累计方差贡献率达到90%时的主成分数量
    # n_components = np.argmax(cumulative_explained_variance >= 0.9) + 1
    # n_components = max(2, n_components)
    # print(f'{n_components} components selected')
    #
    # # n_components = 1
    # svd = TruncatedSVD(n_components=n_components)
    # svd.fit(rt_scores)
    #
    # importance_matrix = np.dot(np.diag(svd.singular_values_), svd.components_)
    # # Compute the L2 norm of each column (Euclidean norm) for "importance"
    # importance = np.linalg.norm(importance_matrix, axis=0)
    # # importance = np.linalg.norm(svd.components_, axis=0)
    #
    # # importance = np.dot(svd.singular_values_, svd.components_)
    # sort_indices = np.argsort(-importance)
    # print("Columns importance: ", list(rt_scores.columns[sort_indices]))
    # print("Importance: ", importance[sort_indices])
    #
    # print("Selected: ", list(rt_scores.columns[sort_indices][: n_components]))
    #
    # print('debug')