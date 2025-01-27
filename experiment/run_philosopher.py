# Download released package from https://github.com/Nesvilab/philosopher
import os
import subprocess
from pathlib import Path

PHILOSOPHER = '/home/nico/software/philosopher/philosopher'

pepxml_files = Path('../data/JY_1_10_25M/no_booster/').glob('*.pep.xml')
output_folder = str(Path('../experiment/JY_1_10_25M/philosopher_no').resolve())

fasta_path = str(Path('../data/JY_1_10_25M/2024-09-03-decoys-contam-Human_EBV_GD1_B95.fasta').resolve())

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

def run_philosopher(pepxml_path):
    file_name = pepxml_path.stem.split('.')[0]
    file_path = pepxml_path.parent
    workspace_name = os.path.basename(output_folder)
    subprocess.run(f'{PHILOSOPHER} workspace --clean', shell=True, cwd=output_folder)
    subprocess.run(f'{PHILOSOPHER} workspace --init', shell=True, cwd=output_folder)
    subprocess.run(f'{PHILOSOPHER} database --annotate {fasta_path}', shell=True, cwd=output_folder)
    # subprocess.run(f'{PHILOSOPHER} peptideprophet --database {fasta_path} --decoyprobs --ppm --accmass --nonparam --expectscore --enzyme nonspecific {pepxml_path}', shell=True, cwd=output_folder)
    # subprocess.run(f'{PHILOSOPHER} filter --pepxml {file_path}/interact-{file_name}.pep.xml', shell=True, cwd=output_folder)
    subprocess.run(f'{PHILOSOPHER} filter --pepxml {file_path}/{file_name}.pep.xml', shell=True, cwd=output_folder)
    subprocess.run(f'{PHILOSOPHER} report --prefix 1', shell=True, cwd=output_folder)
    subprocess.run(f'mv {workspace_name}_ion.tsv {file_name}_ion.tsv', shell=True, cwd=output_folder)
    subprocess.run(f'mv {workspace_name}_psm.tsv {file_name}_psm.tsv', shell=True, cwd=output_folder)
    subprocess.run(f'mv {workspace_name}_peptide.tsv {file_name}_peptide.tsv', shell=True, cwd=output_folder)


if __name__ == '__main__':
    for pepxml_path in pepxml_files:
        run_philosopher(pepxml_path.resolve())
