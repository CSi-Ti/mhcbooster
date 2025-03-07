from pyteomics import mzml
from pathlib import Path
path = Path(__file__)

file = mzml.read('/mnt/d/data/PXD052187/raw/IP0040_11MAI2022_JY_MHC1_HUMAN_S4_2PELLETS_KK_SERIAL_DIL_PT_1_R1_uncalibrated.mzML')

ms2_file = [data for data in file if data['ms level'] == 2]

# import pyopenms
#
# mzml_file = '/mnt/d/data/JY_1_10_25M/timsconvert/JY_Class1_1M_DDA_60min_Slot1-10_1_541.mzML'
# exp = pyopenms.MSExperiment()
# mzml = pyopenms.MzMLFile()
#
# # Load the mzML file into the MSExperiment object
# mzml.load(mzml_file, exp)
# ms2_file = []
# for spectrum in exp.getSpectra():
#     if spectrum.getMSLevel() == 2:
#         ms2_file.append(spectrum)

print('debug')