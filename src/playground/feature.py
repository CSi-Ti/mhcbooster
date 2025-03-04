import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.utils.fdr import calculate_qs
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Lasso

feature_matrix = pd.read_csv('/mnt/d/workspace/mhc-validator-2/experiment/PXD019643/test/160311_DK_AUT01-DN02_Bladder_W6-32_8_DDA_1_400-650mz_msms32_standard_test/160311_DK_AUT01-DN02_Bladder_W6-32_8_DDA_1_400-650mz_msms32_standard.features.tsv', sep='\t')

rt_scores = feature_matrix[[col for col in feature_matrix.columns if 'log_rt_error' in col]]


# Lasso回归进行特征选择
lasso = Lasso(alpha=0.01)
lasso.fit(rt_scores, feature_matrix[['Label']])

# 查看Lasso回归后的系数
print("Lasso回归的系数：", lasso.coef_)

# 根据系数选择最重要的特征
important_features = np.where(lasso.coef_ != 0)[0]
print("选择的特征索引：", important_features)


x_qs = calculate_qs(rt_scores.iloc[:, 0], feature_matrix['Label'], higher_better=False)
y_qs = calculate_qs(rt_scores.iloc[:, 1], feature_matrix['Label'], higher_better=False)


plt.scatter(rt_scores.iloc[:, 0], rt_scores.iloc[:, 1], c=feature_matrix['Label'], s=2)
plt.show()
print()