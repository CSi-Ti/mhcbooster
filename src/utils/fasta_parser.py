
from pyteomics import fasta

def parse_keys(fasta_path):
    fasta_data = fasta.read(fasta_path)
    key_info_map = {}
    for protein in fasta_data:
        description = protein.description
        key = description.split(' ')[0]
        if not key.startswith('VAR|'):
            continue
        if key in key_info_map.keys():
            print('Warning: duplicate key detected in fasta file')
        key_info_map[key] = description
    return key_info_map


if __name__ == '__main__':
    key_info_map = parse_keys('/mnt/d/data/JY_1_10_25M/2025-01-29-decoys-contam-new_modified_fasta_GP2.fasta.fas')
    print('debug')