import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
# import numpy as np
import pandas as pd

# Generate random data for multiple categories
data = pd.read_csv('/mnt/d/workspace/mhc-validator-2/experiment/PXD022194/result_stats.tsv', sep='\t')

percolator = np.ones(len(data))
fragpipe = data['fragpipe'] / data['percolator']
mhcbooster = data['mhcbooster'] / data['percolator']
# Create a swarm plot for the categories
sns.swarmplot(data=pd.DataFrame([percolator, fragpipe, mhcbooster]).T)

# Calculate the mean for each category
# category_means = df.groupby('Category')['Values'].mean()
#
# # Add a horizontal line for each category's mean value
# for category, mean in category_means.items():
#     plt.axhline(mean, color='red', linestyle='--', label=f'Mean ({category}): {mean:.2f}')

# Add labels and title
plt.title("Swarm Plot with Mean Value Lines by Category")
plt.xlabel("Category")
plt.ylabel("Value")

# Show the legend
plt.legend()

# Show the plot
plt.show()
