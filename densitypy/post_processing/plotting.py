import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Palette definitions
palette = [
    [0.00, 0.00, 0.00, 0.50],
    [0.10, 0.00, 0.50, 1.00],
    [0.17, 0.00, 1.00, 1.00],
    [0.25, 0.00, 1.00, 0.50],
    [0.34, 0.50, 1.00, 0.00],
    [0.46, 1.00, 1.00, 0.00],
    [0.61, 1.00, 0.50, 0.00],
    [1.00, 0.50, 0.00, 0.00]
]

# ... you can do similar for other palettes

# plot size
plt.figure(figsize=(53.33, 32))  # size in inches. 3840x2304 pixels = 53.33x32 inches for 72 dpi

# Axes labels
plt.xlabel("Excitation energy (a.u)")
plt.ylabel("Emission energy (a.u)")

# Axes ticks
plt.xticks(np.arange(0, 1.01, 0.04))  # customize as needed
plt.yticks(np.arange(0, 1.01, 0.04))  # customize as needed

# Load and plot data
data = pd.read_csv('/home/ruben/PycharmProjects/DensityPy/Studies/cluttertest/sim/Dipole/DipoleFT_ww',
                   sep='\s+', header=None,
                   names=['X', 'Y', 'X-real-value', 'X-imag-value', 'Y-real-value', 'Y-imag-value', 'Z-real-value',
                          'Z-imag-value'])

# average_ALL = np.sqrt((data[3]**2 + data[4]**2 + data[5]**2 + data[6]**2 + data[7]**2 + data[8]**2) / 3)
#  X    Y  X-real-value  ...  Y-imag-value  Z-real-value  Z-imag-value
average_ww_dim = np.sqrt((data["X-real-value"] ** 2 + data["X-imag-value"] ** 2 + data["Y-real-value"] ** 2 + data[
    "Y-imag-value"] ** 2 + data["Z-real-value"] ** 2 + data["Z-imag-value"] ** 2) / 3)
plt.plot(data["X"], data["Y"], average_ww_dim)  # replace column indices as appropriate

# Save plot
plt.savefig("NMA_DipoleFT_ALL.png")

# ... you can create other plots here as needed

plt.show()
