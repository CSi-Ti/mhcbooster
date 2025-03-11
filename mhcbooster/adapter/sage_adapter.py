import subprocess

import pandas as pd

from pathlib import Path


def split(result_path):
    result_path = Path(result_path)
    result_df = pd.read_csv(result_path, sep='\t')

    if result_df['ln(precursor_ppm)'].isna().any():
        result_df.dropna(subset=['ln(precursor_ppm)'], inplace=True)

    file_names = result_df['FileName'].drop_duplicates().tolist()
    for file_name in file_names:
        run_df = result_df[result_df['FileName'] == file_name]
        run_path = (result_path.parent/file_name).with_suffix('.pin')
        run_df.to_csv(run_path, sep='\t', index=False)
    result_path.rename(result_path.with_suffix('.pin.tmp'))

def get_sage_command(exe_path, param_path, fasta_path, mzml_path, output_path):
    # mhcbooster/third_party/sage/target/release/sage pipeline_setup/sage.json -f /mnt/d/data/JY_1_10_25M/2024-09-03-decoys-contam-Human_EBV_GD1_B95.fasta --write-pin /mnt/d/data/JY_1_10_25M/timsconvert/*.mzML -o experiment/JY_1_10_25M/sage'
    command = f'{exe_path} {param_path} -f {fasta_path} --write-pin {mzml_path} -o {output_path}'
    return command

def run_sage(exe_path, param_path, fasta_path, mzml_path, output_path):
    command = get_sage_command(exe_path, param_path, fasta_path, mzml_path, output_path)
    subprocess.run(command, shell=True)