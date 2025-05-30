# pip install mokapot
import os
from pathlib import Path
from os import system


pin_files = Path('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/sage_test_single/sage').rglob('*.pin')
output_folder = Path('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/sage_test_single/mokapot/')

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for pin in pin_files:
    # FOR SAGE
    # if 'edited' in pin.name:
    #     continue
    # import pandas as pd
    # df = pd.read_csv(pin, sep='\t')
    # df = df.drop(columns=['ln(precursor_ppm)', 'FileName'])
    # pin = str(pin).replace('.pin', '_edited.pin')
    # df.to_csv(pin, sep='\t', index=False)
    # pin = Path(pin)

    file_name = pin.stem
    print(file_name)
    system(f'mokapot --decoy_prefix rev_ -d {output_folder} -r {pin.stem} {pin}')
    # system(f'mokapot -r {pin.stem} {pin}')
