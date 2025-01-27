import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt

# Generate example data (e.g., noisy sinusoidal data)
np.random.seed(42)
x = np.linspace(0, 10, 100)
y = np.sin(x) + np.random.normal(0, 0.1, len(x))  # Noisy sine wave

# Perform LOESS regression using statsmodels' lowess function
frac = 0.2  # Fraction of points used for local regression (controls smoothing)
lowess = sm.nonparametric.lowess(y, x, frac=frac)

# The lowess object contains two columns: x values and smoothed y values
x_smooth = lowess[:, 0]
y_smooth = lowess[:, 1]

# Plot the original data and the LOESS smoothed curve
plt.scatter(x, y, label='Data', color='blue', alpha=0.6)
plt.plot(x_smooth, y_smooth, label='LOESS fit', color='red', linewidth=2)
plt.legend()
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f'LOESS Regression (frac={frac})')
plt.show()