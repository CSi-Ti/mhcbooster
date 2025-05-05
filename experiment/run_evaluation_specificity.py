import numpy as np
import pandas as pd

from pathlib import Path
import matplotlib.pyplot as plt
from mhcbooster.utils.peptide import remove_modifications, remove_charge, remove_previous_and_next_aa

def get_pin_r_map(pin_folder):
    def get_r_from_proteins_list(proteins_list):
        n_psm_yeast = 0
        for proteins in proteins_list:
            prots = proteins.split(';')
            yeast, human = False, False
            for prot in prots:
                if prot.startswith('rev_'):
                    continue
                if 'YEAST' in prot:
                    yeast = True
                if 'HUMAN' in prot:
                    human = True
            if yeast and not human:
                n_psm_yeast += 1
        r = n_psm_yeast / (len(proteins_list) - n_psm_yeast)
        return r

    pin_r_map = {}
    pin_paths = Path(pin_folder).rglob('*.pin')
    for pin in pin_paths:
        file_name = pin.stem
        pin_df = pd.read_csv(pin, sep='\t')
        pin_df = pin_df[pin_df['Label'] == 1]
        pin_df = pin_df[['Peptide', 'Proteins']]
        pin_df['peptide'] = remove_charge(remove_previous_and_next_aa(pin_df['Peptide'].to_numpy()))
        pin_df['sequence'] = remove_modifications(pin_df['peptide'].to_numpy())
        peptide_df = pin_df.groupby('peptide', as_index=False).agg({'Proteins': 'first'})
        sequence_df = pin_df.groupby('sequence', as_index=False).agg({'Proteins': 'first'})

        r_psm = get_r_from_proteins_list(pin_df['Proteins'].to_numpy())
        r_peptide = get_r_from_proteins_list(peptide_df['Proteins'].to_numpy())
        r_sequence = get_r_from_proteins_list(sequence_df['Proteins'].to_numpy())
        pin_r_map[file_name] = [r_psm, r_peptide, r_sequence]
    return pin_r_map

def generate_fdp_data(data_df, r):
    n_human, n_yeast = 0, 0
    fdp_data = []
    for index, row in data_df.iterrows():
        proteins = row['proteinIds']
        qvalue = row['q-value']
        prots = proteins.split(';')
        if len(prots) == 1:
            prots = proteins.split(',')
        yeast, human = False, False
        for prot in prots:
            if prot.startswith('rev_'):
                continue
            if 'YEAST' in prot:
                yeast = True
            if 'HUMAN' in prot:
                human = True
        if yeast and not human:
            n_yeast += 1
        else:
            n_human += 1
        lower_bound = n_yeast / (n_yeast + n_human)
        upper_bound = n_yeast * (1 + 1 / r) / (n_yeast + n_human)
        fdp_data.append([qvalue, lower_bound, upper_bound])
    return fdp_data

def draw_fdp_distribution(fdp_data, file_name):
    fdp_data = np.array(fdp_data)
    plt.plot(fdp_data[:, 0], fdp_data[:, 0], c='g')
    plt.plot(fdp_data[:, 0], fdp_data[:, 1], c='r', label='FDP lower bound')
    plt.plot(fdp_data[:, 0], fdp_data[:, 2], c='b', label='FDP upper bound')
    plt.xlim(0, 0.1)
    plt.ylim(0, 0.1)
    plt.legend()
    plt.title(file_name)
    plt.xlabel('FDR threshold')
    plt.ylabel('Estimated FDP')
    plt.show()

def eval_pout(result_folder, pin_r_map, tool_name):
    def generate_modified_pout(path):
        # Compatibility for multi-prot
        with open(path, 'r') as file:
            content = file.read().replace(';\t', ';')
        with open(str(path) + '.m', 'w') as file:
            file.write(content)

    psm_paths = Path(result_folder).rglob('*_psm_target.pout')
    for psm_path in psm_paths:
        file_name = str(psm_path.name).replace('_psm_target.pout', '')
        generate_modified_pout(str(psm_path))
        psm_df = pd.read_csv(str(psm_path) + '.m', sep='\t')
        psm_df = psm_df[['proteinIds', 'q-value']]
        psm_r = pin_r_map[file_name][0]
        psm_fdp_data = generate_fdp_data(psm_df, psm_r)
        draw_fdp_distribution(psm_fdp_data, file_name + f'_{tool_name}_psm')

    pep_paths = Path(result_folder).rglob('*_pep_target.pout')
    for pep_path in pep_paths:
        file_name = str(pep_path.name).replace('_pep_target.pout', '')
        generate_modified_pout(str(pep_path))
        pep_df = pd.read_csv(str(pep_path) + '.m', sep='\t')

        pep_df['pep'] = remove_charge(remove_previous_and_next_aa(pep_df['peptide'].to_numpy()))
        pep_df['seq'] = remove_modifications(pep_df['pep'].to_numpy())

        peptide_df = pep_df.groupby('pep', as_index=False).agg({'proteinIds': 'first', 'q-value': 'first'})
        peptide_df = peptide_df.sort_values(by='q-value')
        pep_r = pin_r_map[file_name][1]
        pep_fdp_data = generate_fdp_data(peptide_df, pep_r)
        draw_fdp_distribution(pep_fdp_data, file_name + f'_{tool_name}_pep')

        sequence_df = pep_df.groupby('seq', as_index=False).agg({'proteinIds': 'first', 'q-value': 'first'})
        sequence_df = sequence_df.sort_values(by='q-value')
        seq_r = pin_r_map[file_name][2]
        seq_fdp_data = generate_fdp_data(sequence_df, seq_r)
        draw_fdp_distribution(seq_fdp_data, file_name + f'_{tool_name}_seq')

def eval_tsv(result_folder, pin_r_map, tool_name):
    psm_paths = Path(result_folder).rglob('psm.tsv')
    for psm_path in psm_paths:
        file_name = str(psm_path.parent.name).replace('_MHCBooster', '')
        psm_df = pd.read_csv(psm_path, sep='\t')
        psm_df = psm_df[psm_df['label'] == 'Target']
        # psm_df = psm_df[psm_df['psm_qvalue'] <= 0.01]
        psm_df['proteinIds'] = psm_df['protein'].astype(str) + ',' + psm_df['mapped_protein'].astype(str)
        psm_df['q-value'] = psm_df['psm_qvalue']
        psm_df = psm_df.sort_values(by='q-value')
        psm_r = pin_r_map[file_name][0]
        psm_fdp_data = generate_fdp_data(psm_df, psm_r)
        draw_fdp_distribution(psm_fdp_data, file_name + f'_{tool_name}_psm')

    pep_paths = Path(result_folder).rglob('peptide.tsv')
    for pep_path in pep_paths:
        file_name = str(pep_path.parent.name).replace('_MHCBooster', '')
        pep_df = pd.read_csv(pep_path, sep='\t')
        pep_df = pep_df[pep_df['label'] == 'Target']
        # pep_df = pep_df[pep_df['pep_qvalue'] <= 0.01]
        pep_df['proteinIds'] = pep_df['protein'].astype(str) + ',' + pep_df['mapped_protein'].astype(str)
        pep_df['q-value'] = pep_df['pep_qvalue']
        pep_df = pep_df.sort_values(by='q-value')
        pep_r = pin_r_map[file_name][1]
        pep_fdp_data = generate_fdp_data(pep_df, pep_r)
        draw_fdp_distribution(pep_fdp_data, file_name + f'_{tool_name}_pep')

        sequence_df = pep_df.groupby('sequence', as_index=False).agg({'proteinIds': 'first', 'q-value': 'first'})
        sequence_df = sequence_df.sort_values(by='q-value')
        seq_r = pin_r_map[file_name][2]
        seq_fdp_data = generate_fdp_data(sequence_df, seq_r)
        draw_fdp_distribution(seq_fdp_data, file_name + f'_{tool_name}_seq')

def count_yeast_matches(result_folder, suffix):
    if '.pout' in suffix:
        q_col = 'q-value'
        label_col = None
        label_target = None
        prot_cols = ['proteinIds']
    else:
        q_col = 'pep_qvalue'
        label_col = 'label'
        label_target = 'Target'
        prot_cols = ['protein', 'mapped_protein']

    result_files = Path(result_folder).rglob('*' + suffix)
    for result_file in result_files:
        if 'combined' in result_file.name:
            continue
        result_df = pd.read_csv(result_file, sep='\t')
        if label_col is not None:
            result_df= result_df[result_df[label_col] == label_target]
        result_df = result_df[result_df[q_col] <= 0.01]
        result_df = result_df.reset_index(drop=True)
        peptides = remove_charge(remove_previous_and_next_aa(result_df['peptide'].to_numpy()))
        sequences = remove_modifications(peptides)

        result_df = result_df.iloc[[i for i in range(len(result_df)) if (8 <= len(sequences[i]) <= 15)]]
        result_df = result_df.reset_index(drop=True)


        peptide_num = len(np.unique(peptides))
        sequence_num = len(np.unique(sequences))
        yeast_peps = []
        for i in range(len(result_df)):
            prot_row = result_df.loc[i, prot_cols]
            prots = ';'.join([p for p in prot_row if pd.notna(p)])
            if 'YEAST' in prots and 'HUMAN' not in prots:
                yeast_peps.append(peptides[i])
        yeast_num = len(np.unique(yeast_peps))
        if '.pout' in suffix:
            print(result_file.name.replace(suffix, ''))
        else:
            print(result_file.parent.name.replace('_MHCBooster', ''))
        print(f'{peptide_num}\t{sequence_num}\t{yeast_num}\t{yeast_num/peptide_num}')


if __name__ == '__main__':
    pin_folder = Path('/mnt/d/data/JY_1_10_25M/fragpipe/Search_human_yeast')
    result_folder = Path('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/human_yeast')
    percolator_result_folder = result_folder / 'percolator'
    fragpipe_result_folder = result_folder / 'fragpipe'
    mhcboost_result_folder = result_folder / 'mhcbooster'

    # count_yeast_matches(percolator_result_folder, '_pep_target.pout')
    # count_yeast_matches(fragpipe_result_folder, '_edited_pep_target.pout')
    # count_yeast_matches(mhcboost_result_folder, 'peptide.tsv')
    count_yeast_matches(result_folder / 'test_nonorm_specificity', 'peptide.tsv')
    print('--------')
    # count_yeast_matches(result_folder / 'test_norm', 'peptide.tsv')
    # count_yeast_matches(result_folder / 'test_norm_re', 'peptide.tsv')
    # count_yeast_matches(result_folder / 'test_no_norm', 'peptide.tsv')
    # count_yeast_matches(result_folder / 'mhcbooster_test_focal_0alpha_oldk_3', 'peptide.tsv')

    # pin_r_map = get_pin_r_map(pin_folder)
    # eval_pout(percolator_result_folder, pin_r_map, 'perc')
    # eval_pout(fragpipe_result_folder, pin_r_map, 'frag')
    # eval_tsv(result_folder / 'test_nonorm_specificity', pin_r_map, 'mhcb')
    # print(pin_r_map)
