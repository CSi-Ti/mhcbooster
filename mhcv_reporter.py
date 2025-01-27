import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def read_result_file(path, feature_map):
    mhcv_df = pd.read_csv(path, sep='\t')
    mhcv_df['ptm'] = mhcv_df['Peptide'].apply(lambda pep: pep[2:-3])
    target_df = mhcv_df[mhcv_df['mhcv_label'] == 1]
    psm_df = target_df[target_df['mhcv_q-value'] <= 0.01]
    pep_df = target_df[target_df['mhcv_pep-level_q-value'] <= 0.01]
    pep_df = pep_df.drop_duplicates(subset=['ptm', 'Proteins'])


    peptides = pep_df[['ptm', 'mhcv_peptide', 'Proteins']].copy()
    peptides = peptides.groupby('ptm', as_index=False).agg({'mhcv_peptide': lambda x: x.iloc[0],
                                                            'Proteins': lambda x: ''.join(x)})

    # Append organism
    def get_org_from_prot(row):
        orgs = set()
        for prot in row['Proteins'].split(';'):
            if prot == '' or 'rev_' in prot:
                continue
            prot_split = prot.split('_')
            if len(prot_split) > 1:
                orgs.add(prot.split('_')[-1])
        return ';'.join(orgs)
    peptides['Organism'] = peptides.apply(get_org_from_prot, axis=1)

    # Append binding status
    netmhcpan_path = str(path).replace('MhcValidator_annotated', 'NetMHCpan_Predictions')
    if os.path.exists(netmhcpan_path):
        netmhcpan_df = pd.read_csv(netmhcpan_path, sep='\t')
        netmhcpan_df = netmhcpan_df[['Peptide', 'Binder']]
        netmhcpan_df = netmhcpan_df.drop_duplicates(['Peptide'])
        netmhcpan_df = netmhcpan_df.rename(columns={'Peptide': 'mhcv_peptide'})
        peptides = peptides.merge(netmhcpan_df, how='left', on='mhcv_peptide')
    else:
        peptides['Binder'] = ''

    # Append features (RT, MS2, CCS, etc.)
    if len(feature_map) > 0:
        pep_feature_df = (pep_df[['ptm'] + list(feature_map.keys())].groupby(['ptm'], as_index=False)
                          .agg(feature_map))
        peptides = peptides.merge(pep_feature_df, how='left', on='ptm')

    # Append matched psm number
    pep_stat = target_df['ptm'].value_counts().reset_index()
    peptides = peptides.merge(pep_stat, how='left', on='ptm')
    file_name = str(path.name).split('.')[0]

    # Set column names
    peptides.columns = ['Peptide', 'Sequence', 'Proteins', 'Organism', 'Binder'] + list(feature_map.keys()) + [file_name + ' Spectral Count']

    return len(psm_df), peptides

def combine_result_files(root_folder, feature_map):
    mhcv_result_paths = root_folder.rglob('*.MhcValidator_annotated.tsv')

    ident_numbers = []
    result_table = pd.DataFrame(columns=['Peptide', 'Sequence', 'Proteins', 'Organism', 'Binder'])
    for path in mhcv_result_paths:
        file_name = str(path.name).split('.')[0]
        n_psm, peptides = read_result_file(path, feature_map)

        n_strong = np.sum(peptides['Binder'] == 'Strong')
        n_weak = np.sum(peptides['Binder'] == 'Weak')
        n_none = np.sum(peptides['Binder'] == 'Non-binder')
        ident_numbers.append([file_name, n_psm, len(peptides), n_strong, n_weak, n_none])

        result_table = result_table.merge(peptides, on=['Peptide', 'Sequence', 'Proteins', 'Organism', 'Binder'], how='outer',
                                          suffixes=('_df1', '_df2'))
        if len(feature_map) != 0 and list(feature_map.keys())[0] not in result_table.columns:
            for feature in feature_map.keys():
                if feature_map[feature] == 'max':
                    result_table[feature] = result_table[[f'{feature}_df1', f'{feature}_df2']].max(axis=1)
                else:
                    result_table[feature] = result_table[[f'{feature}_df1', f'{feature}_df2']].min(axis=1)
                result_table = result_table.drop(columns=[f'{feature}_df1', f'{feature}_df2'])
    ident_numbers = pd.DataFrame(ident_numbers,
                                 columns=['File name', 'PSMs', 'Peptides', 'Strong binders', 'Weak binders',
                                          'None binders'])
    ident_numbers.to_csv(root_folder / 'combined_result_stat.tsv', sep='\t', index=False)
    result_table.fillna(0, inplace=True)
    if len(feature_map) > 0:
        result_table = result_table[result_table.columns[:5].append(result_table.columns[-len(feature_map):]).append(
            result_table.columns[5:-len(feature_map)])]
    result_table[result_table.columns[5 + len(feature_map):]] = result_table[
        result_table.columns[5 + len(feature_map):]].astype(int)

    agg_dict = {result_table.columns[i]: 'first' for i in [1] + list(range(4, 5+len(feature_map)))}
    agg_dict.update({
        result_table.columns[2]: lambda x: ''.join(x),
        result_table.columns[3]: lambda x: ''.join(x),
    })
    agg_dict.update({
        result_table.columns[i]: 'sum' for i in range(5+len(feature_map), len(result_table.columns))
    })
    result_table = result_table.groupby('Peptide', as_index=False).agg(agg_dict)

    def protein_filtering(x):
        x_combine = set()
        for split in x.split(';'):
            if split == '' or split.startswith('rev_'):
                continue
            x_combine.add(split.replace('@@', ''))
        return ';'.join(sorted(x_combine))

    result_table['Proteins'] = result_table['Proteins'].apply(protein_filtering)

    result_table.to_csv(root_folder / 'combined_result.tsv', sep='\t', index=False)
    print('Combined result generated.')
    return result_table

def draw_feature_distribution(root_folder, feature_map, organism):
    result_table = pd.read_csv(root_folder / 'combined_result.tsv', sep='\t')
    result_table_binder = result_table[(result_table['Binder'] == 'Strong') | (result_table['Binder'] == 'Weak')]

    for feature in feature_map.keys():
        counts, bins = np.histogram(result_table_binder[feature], bins=50)
        normalized_counts = counts / counts.sum()
        plt.bar(bins[:-1], normalized_counts, width=np.diff(bins), alpha=0.7, color='blue', label='All binders')

        organism_binder = result_table_binder[result_table_binder['Organism'].str.contains(organism)]
        counts, bins = np.histogram(organism_binder[feature], bins=bins)
        normalized_counts = counts / counts.sum()
        plt.bar(bins[:-1], normalized_counts, width=np.diff(bins), alpha=0.7, color='orange', label=f'{organism} binders')

        plt.title(feature)
        plt.legend()
        plt.xlabel(feature)
        plt.ylabel('Percentage')
        # plt.show()
        plt.savefig(root_folder / f'{organism}_distribution_{feature}.png')
        plt.close()


if __name__ == '__main__':
    root_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/JPST002044/msb')
    feature_map = {
        # 'unweighted_spectral_entropy': 'max',
        # 'delta_RT_loess': 'min',
        # 'im2deep_delta_ccs': 'min',
        # 'autort_delta_rt': 'min'
    }
    combine_result_files(root_folder, feature_map)
    # draw_feature_distribution(root_folder, feature_map, 'EBVB9')
    # draw_feature_distribution(root_folder, feature_map, 'HHV8P')
    # draw_feature_distribution(root_folder, feature_map, 'MLVFF')
