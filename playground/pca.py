import numpy as np
import pandas as pd

from mhcvalidator.fdr import calculate_qs
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

feature_matrix = pd.read_csv('/mnt/d/workspace/mhc-validator-2/experiment/PXD019643/test/160311_DK_AUT01-DN02_Bladder_W6-32_8_DDA_1_400-650mz_msms32_standard_test/160311_DK_AUT01-DN02_Bladder_W6-32_8_DDA_1_400-650mz_msms32_standard.features.tsv', sep='\t')

rt_scores = feature_matrix[[col for col in feature_matrix.columns if 'log_rt_error' in col]]
# 将数据转换为 DataFrame
df = pd.DataFrame(rt_scores)

# 标准化数据
scaler = StandardScaler()
scaled_data = scaler.fit_transform(df)

# 应用 PCA
pca = PCA(n_components=4)  # 假设保留所有4个主成分
pca.fit(scaled_data)

# 输出主成分贡献的方差比例
print("各主成分的方差解释比例:", pca.explained_variance_ratio_)

# 输出每个特征在主成分中的权重
print("\n主成分的每个特征的权重 (components_):")
print(pd.DataFrame(pca.components_, columns=df.columns))

# 根据每个主成分的特征权重，选择最重要的特征
# 比较每个特征在所有主成分中的贡献
feature_importance = np.abs(pca.components_).sum(axis=0)

# 创建一个 DataFrame 来查看每个特征的总贡献度
feature_importance_df = pd.DataFrame({'Feature': df.columns, 'Importance': feature_importance})
feature_importance_df = feature_importance_df.sort_values(by='Importance', ascending=False)

print("\n根据主成分贡献度排序的特征重要性:")
print(feature_importance_df)