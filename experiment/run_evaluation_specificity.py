import numpy as np
import pandas as pd

from pathlib import Path
from mhcbooster.utils.peptide import remove_modifications, remove_charge, remove_previous_and_next_aa

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
        sequences = remove_modifications(remove_charge(remove_previous_and_next_aa(result_df['peptide'].to_numpy())))
        result_df = result_df.iloc[[i for i in range(len(result_df)) if (8 <= len(sequences[i]) <= 15)]]
        result_df = result_df.reset_index(drop=True)

        peptide_num = len(result_df)
        yeast_num = 0
        for i in range(len(result_df)):
            prot_row = result_df.loc[i, prot_cols]
            prots = ';'.join([p for p in prot_row if pd.notna(p)])
            if 'YEAST' in prots and 'HUMAN' not in prots:
                yeast_num += 1
        if '.pout' in suffix:
            print(result_file.name.replace(suffix, ''))
        else:
            print(result_file.parent.name.replace('_MHCBooster', ''))
        print(f'{peptide_num}\t{len(np.unique(sequences))}\t{yeast_num}\t{yeast_num/peptide_num}')


if __name__ == '__main__':
    result_folder = Path('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/human_yeast')
    percolator_result_folder = result_folder / 'percolator'
    fragpipe_result_folder = result_folder / 'fragpipe'
    mhcboost_result_folder = result_folder / 'mhcbooster'

    # count_yeast_matches(percolator_result_folder, '_pep_target.pout')
    count_yeast_matches(fragpipe_result_folder, '_edited_pep_target.pout')
    # count_yeast_matches(mhcboost_result_folder, 'peptide.tsv')
    count_yeast_matches(result_folder / 'test_norm', 'peptide.tsv')
    count_yeast_matches(result_folder / 'test_norm_re', 'peptide.tsv')
    count_yeast_matches(result_folder / 'test_no_norm', 'peptide.tsv')
    count_yeast_matches(result_folder / 'mhcbooster_test_focal_0alpha_oldk_3', 'peptide.tsv')
