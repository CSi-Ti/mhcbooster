import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyteomics.parser import length
from src.utils.peptide import replace_uncommon_aas, remove_charge, remove_previous_and_next_aa, remove_modifications
from src.predictors.netmhcpan_helper import NetMHCpanHelper

N_THREADS = os.cpu_count() // 2


class Evaluation:
    def __init__(self, dataset_name, min_len, max_len):
        self.result_folder = Path(__file__).parent/dataset_name
        self.min_len = min_len
        self.max_len = max_len

    def base_eval_pep(self, folder, pep_file_suffix, sep, pep_col, qvalue_col, label_col=None, target_label=None, mod_col=None):
        paths = list(Path(folder).rglob(pattern='*' + pep_file_suffix))
        paths.sort()

        result_stat = pd.DataFrame(columns=['File name', 'Peptides', 'Peptides (filtered)'])
        for i, path in enumerate(paths):
            filename = path.stem.replace('_edited_pep_target.pout', '').replace('_pep_target.pout', '').replace('.MhcValidator_annotated', '').replace('_result', '')
            if length(list(open(path))) == 0:
                result_stat.loc[i] = [filename, 0, 0]
                continue

            pep_df = pd.read_csv(path, sep=sep)
            if label_col and target_label:
                pep_df = pep_df[pep_df[label_col] == target_label]
            if qvalue_col:
                pep_df = pep_df[pep_df[qvalue_col] <= 0.01]

            peptides = pep_df[pep_col].to_numpy()
            peptides = remove_previous_and_next_aa(peptides)
            peptides = remove_charge(peptides)
            if mod_col:
                pep_mod_keys = np.char.add(peptides.astype(str), pep_df[mod_col].fillna('').to_numpy().astype(str))
                unique_indices = np.unique(pep_mod_keys, return_index=True)[1]
                peptides = peptides[unique_indices]
            else:
                peptides = np.unique(peptides)

            peptides = remove_modifications(peptides)
            pep_filtered = np.array(peptides)
            if len(pep_filtered) > 0:
                pep_filtered = pep_filtered[
                    (np.char.str_len(pep_filtered) >= self.min_len) * (np.char.str_len(pep_filtered) <= self.max_len)]
            result_stat.loc[i] = [filename, len(pep_filtered), len(np.unique(pep_filtered))]
        return result_stat

    def eval_mokapot(self):
        mokapot_folder = self.result_folder/'mokapot'
        if mokapot_folder.exists():
            result_df = self.base_eval_pep(mokapot_folder, pep_file_suffix='.mokapot.peptides.txt', sep='\t', pep_col='Peptide',
                                      qvalue_col='mokapot q-value', label_col='Label', target_label=True, mod_col=None)
            result_df.columns = ['File name', 'mokapot_pep', 'mokapot_seq']
            return result_df

    def eval_percolator(self):
        percolator_folder = self.result_folder/'percolator'
        if percolator_folder.exists():
            paths = Path(percolator_folder).rglob(pattern='*' + '_pep_target.pout')
            for i, path in enumerate(paths):
                # Compatibility for multi-prot
                with open(path, 'r') as file:
                    content = file.read().replace(';\t', ';')
                with open(str(path) + '.m', 'w') as file:
                    file.write(content)
            result_df = self.base_eval_pep(percolator_folder, pep_file_suffix='_pep_target.pout.m', sep='\t', pep_col='peptide',
                                      qvalue_col='q-value', label_col=None, mod_col=None)
            result_df.columns = ['File name', 'percolator_pep', 'percolator_seq']
            return result_df

    def eval_fragpipe(self):
        fragpipe_folder = self.result_folder/'fragpipe'
        if fragpipe_folder.exists():
            paths = Path(fragpipe_folder).rglob(pattern='*' + '_pep_target.pout')
            for i, path in enumerate(paths):
                # Compatibility for multi-prot
                with open(path, 'r') as file:
                    content = file.read().replace(';\t', ';')
                with open(str(path) + '.m', 'w') as file:
                    file.write(content)
            result_df = self.base_eval_pep(fragpipe_folder, pep_file_suffix='_pep_target.pout.m', sep='\t', pep_col='peptide',
                                      qvalue_col='q-value', label_col=None, mod_col=None)
            result_df.columns = ['File name', 'fragpipe_pep', 'fragpipe_seq']
            return result_df

    def eval_philosopher(self):
        philosopher_folder = self.result_folder/'philosopher'
        if philosopher_folder.exists():
            result_df = self.base_eval_pep(philosopher_folder, pep_file_suffix='_peptide.tsv', sep='\t', pep_col='Peptide',
                                      qvalue_col=None, label_col=None, mod_col='Assigned Modifications')
            result_df.columns = ['File name', 'philosopher_pep', 'philosopher_seq']
            return result_df

    def eval_ms2rescore(self):
        ms2rescore_folder = self.result_folder/'ms2rescore'
        if ms2rescore_folder.exists():
            result_df = self.base_eval_pep(ms2rescore_folder, pep_file_suffix='_result.csv', sep=',', pep_col='peptidoform',
                                           qvalue_col='peptide_qvalue', label_col='is_decoy', target_label=False, mod_col=None)
            result_df.columns = ['File name', 'ms2rescore_pep', 'ms2rescore_seq']
            return result_df

    def eval_mhcbooster(self):
        mhcbooster_folder = self.result_folder/'mhcbooster'
        if mhcbooster_folder.exists():
            result_df = self.base_eval_pep(mhcbooster_folder, pep_file_suffix='.MhcValidator_annotated.tsv', sep='\t', pep_col='Peptide',
                                      qvalue_col='mhcv_pep-level_q-value', label_col='mhcv_label', target_label=1, mod_col=None)
            result_df.columns = ['File name', 'mhcbooster_pep', 'mhcbooster_seq']
            return result_df

    def run(self):
        percolator_df = self.eval_percolator()
        philosopher_df = self.eval_philosopher()
        mokapot_df = self.eval_mokapot()
        ms2rescore_df = self.eval_ms2rescore()
        fragpipe_df = self.eval_fragpipe()
        mhcbooster_df = self.eval_mhcbooster()
        result_dfs = [percolator_df, philosopher_df, mokapot_df, ms2rescore_df, fragpipe_df, mhcbooster_df]

        result_df = pd.DataFrame()
        for df in result_dfs:
            if df is not None:
                result_df = pd.merge(result_df, df, on='File name', how='outer') if not result_df.empty else df
        result_df.to_csv(os.path.join(self.result_folder, 'result_stats.tsv'), sep='\t', index=False)
        print(result_df)

if __name__ == '__main__':

    evaluation = Evaluation('JY_500M', min_len=8, max_len=15)
    evaluation.run()
    evaluation = Evaluation('JY_Fractionation_Replicate_1', min_len=8, max_len=15)
    evaluation.run()
    evaluation = Evaluation('RA_Fractionation_Replicate_1', min_len=8, max_len=15)
    evaluation.run()

    MIN_LENGTH = 9
    MAX_LENGTH = 25
