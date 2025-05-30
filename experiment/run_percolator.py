# conda install bioconda::percolator
import os
from pathlib import Path
from os import system

pin_files = Path('/mnt/d/data/RA_Fractionation/Search_0523').rglob('*.pin')
output_folder = Path('/mnt/d/workspace/mhc-booster/experiment/RA_Fractionation/percolator')

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for pin in pin_files:
    file_name = pin.stem
    if 'edited' in file_name:
        continue
    print(file_name)
    psm_target_out = output_folder / (file_name + '_psm_target.pout')
    psm_decoy_out = output_folder / (file_name + '_psm_decoy.pout')
    pep_target_out = output_folder / (file_name + '_pep_target.pout')
    pep_decoy_out = output_folder / (file_name + '_pep_decoy.pout')
    # system(f'percolator -r {pep_target_out} -B {pep_decoy_out} -m {psm_target_out} -M {psm_decoy_out} {pin}')
    system(f'percolator -r {pep_target_out} -m {psm_target_out} {pin}')
