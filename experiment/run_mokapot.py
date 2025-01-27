# pip install mokapot
import os
from pathlib import Path
from os import system

# pin_files = Path('../data/JY_1_10_25M/no_booster/').glob('*.pin')
pin_files = Path('/mnt/d/data/JY_1_10_25M/fragpipe/Search_MSBooster_FDR1').rglob('*edited_im2deep.pin')
output_folder = '../experiment/JY_1_10_25M/mokapot_ccs'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for pin in pin_files:
    file_name = pin.stem
    print(file_name)
    system(f'mokapot --decoy_prefix rev_ -d {output_folder} -r {pin.stem} {pin}')
    # system(f'mokapot -r {pin.stem} {pin}')
