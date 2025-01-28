import random

import numpy as np
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d



# Generate some example data
np.random.seed(0)
x = -np.linspace(0, 10, 100)
y = np.sin(x) + 0.3 * np.random.randn(100)  # Sine with noise

idx = random.shuffle

# Perform LOESS regression
smoothed = lowess(y, x, frac=0.2)  # frac controls the smoothing span

# Extract the smoothed values
x_smooth = smoothed[:, 0]
y_smooth = smoothed[:, 1]

# Plot the original data and the smoothed curve
plt.scatter(x, y, label='Original Data')
plt.plot(x_smooth, y_smooth, color='red', label='LOESS Smoothing')
plt.legend()
plt.show()


# Create a piecewise linear function (linear spline) from smoothed data
linear_spline = interp1d(x_smooth, y_smooth, kind='linear', fill_value="extrapolate")

# Generate new x values and use the spline to get corresponding y values
x_new = np.linspace(0, 10, 200)
y_new = linear_spline(x_new)

# Plot the linear spline along with the original and smoothed data
plt.scatter(x, y, label='Original Data')
plt.plot(x_smooth, y_smooth, color='red', label='LOESS Smoothing')
plt.plot(x_new, y_new, color='green', label='Linear Spline')
plt.legend()
plt.show()

from scipy.interpolate import CubicSpline

# Example data (noisy sine wave)
x = np.linspace(0, 10, 10)   # 10 data points
y = np.sin(x) + 0.1 * np.random.randn(10)  # Sine curve with noise

# Fit a cubic spline to the data with natural boundary conditions (default)
cs = CubicSpline(x, y)

# Create new x values to evaluate the spline at
x_new = np.linspace(0, 10, 1000)
y_new = cs(x_new)  # Interpolated values at new x points

# Plot the original data and the cubic spline
plt.scatter(x, y, label="Data Points", color='red')
plt.plot(x_new, y_new, label="Cubic Spline", color='blue')
plt.legend()
plt.show()
