
import matplotlib.pyplot as plt
import numpy as np
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mass

# Retrieve the spectrum by its USI.
peptide = "AAAAAAAAAAAAA"
peptide_mass = mass.calculate_mass(peptide)
charge_state = 3
peptide_mz = (peptide_mass + (charge_state * 1.007276)) / charge_state
mzs_top = np.array([214.11862, 285.15573, 356.19284, 427.22995, 498.26706, 569.3042, 640.3413, 711.3784, 782.4155, 853.45264, 303.1663, 374.2034, 445.2405, 516.27765, 587.31476, 658.35187, 729.389])
ints_top = np.array([0.058000322, 0.11392251, 0.13280249, 0.14267457, 0.16630515, 0.13712373, 0.13238974, 0.055706777, 0.02375361, 0.004642043, 0.0011714724, 0.0053157704, 0.008672275, 0.0057611926, 0.007857157, 0.0030173366, 0.00088385516])
ints_top = ints_top[np.argsort(mzs_top)]
mzs_top = np.sort(mzs_top)
spectrum_top = sus.MsmsSpectrum(peptide, peptide_mz, charge_state, mzs_top, ints_top, 1)
spectrum_top = spectrum_top.annotate_proforma(peptide, 20, "ppm")

mzs_bottom = np.array([72.04439, 90.054955, 143.0815, 161.09207, 214.11862, 232.12918, 285.15573, 303.1663, 356.19284, 374.2034, 427.22995, 445.2405, 498.26706, 516.27765, 569.3042, 587.31476, 640.3413, 658.35187, 711.3784, 729.389, 782.4155, 800.4261, 853.45264, 871.4632])
ints_bottom = np.array([0.000585687, 0.011478545, 0.08412633, 0.074349344, 0.098416984, 0.0050035985, 0.07764679, 0.0024212322, 0.093568265, 0.002461908, 0.121157214, 0.002367523, 0.12113316, 0.0017787149, 0.14215595, 0.0011818329, 0.10583924, 0.0013559341, 0.03691124, 0.0017238454, 0.011039282, 0.0001858432, 0.0026443356, 0.00046720088])
spectrum_bottom = sus.MsmsSpectrum(peptide, peptide_mz, charge_state, mzs_bottom, ints_bottom, 1)
spectrum_bottom = spectrum_bottom.annotate_proforma(peptide, 20, "ppm")

# Plot the spectrum.
fig, ax = plt.subplots(figsize=(12, 6))
sup.mirror(spectrum_top, spectrum_bottom, ax=ax)
# sup.spectrum(spectrum, grid=False, ax=ax)
ax.set_title(peptide, fontdict={"fontsize": "xx-large"})
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
plt.savefig("proforma_ex3.png", bbox_inches="tight", dpi=300, transparent=False)
plt.close()