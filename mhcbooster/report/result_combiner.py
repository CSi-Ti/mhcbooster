
import pandas as pd
from pathlib import Path

class ResultCombiner:

    def __init__(self, result_folder, psm_fdr=0.01, pep_fdr=0.01, remove_contaminant=False):

        self.result_folder = result_folder
        self.psm_fdr = psm_fdr
        self.pep_fdr = pep_fdr
        self.decoy_prefix = 'rev_'

    def run(self):
        print('Combining results and generating reports...')
        self.psm_result_map = self.read_psm_results()
        self.pep_result_map = self.psm_to_pep(ignore_mod=False)
        self.seq_result_map = self.psm_to_pep(ignore_mod=True)
        self.combined_pep = self.combine_pep(self.pep_result_map, ignore_mod=False)
        self.combined_seq = self.combine_pep(self.seq_result_map, ignore_mod=True)

        self.combined_pep.to_csv(Path(self.result_folder, 'combined_peptides.tsv'), sep='\t', index=False)
        self.combined_seq.to_csv(Path(self.result_folder, 'combined_sequences.tsv'), sep='\t', index=False)
        print(f'Reports saved to {self.result_folder}')


    def read_psm_results(self):
        psm_result_map = {}
        result_paths = Path(self.result_folder).rglob('*.MhcValidator_annotated.tsv')
        for path in result_paths:
            result_df = pd.read_csv(path, sep='\t')
            formatted_df = pd.DataFrame()
            formatted_df['peptide'] = result_df['Peptide'].apply(lambda pep: pep[2:-3])
            formatted_df['sequence'] = result_df['mhcv_peptide']
            formatted_df['prev_AA'] = result_df['Peptide'].apply(lambda pep: pep[0])
            formatted_df['next_AA'] = result_df['Peptide'].apply(lambda pep: pep[-1])
            formatted_df['label'] = result_df['mhcv_label'].map({1: 'Target', 0: 'Decoy'})
            formatted_df['seq_len'] = formatted_df['sequence'].str.len().astype(str)
            formatted_df['charge'] = result_df['Peptide'].apply(lambda pep: str(pep[-3]))
            formatted_df['score'] = result_df['mhcv_prob']
            formatted_df['psm_qvalue'] = result_df['mhcv_q-value']
            formatted_df['pep_qvalue'] = result_df['mhcv_pep-level_q-value']
            formatted_df['protein'] = result_df['Proteins'].str.replace('@', '', regex=False)
            psm_result_map[path.name.split('.')[0]] = formatted_df
        return psm_result_map

    def psm_to_pep(self, ignore_mod=False):
        pep_result_map = {}
        for file in self.psm_result_map.keys():
            psm_df = self.psm_result_map[file]
            psm_df = psm_df[psm_df['label'] == 'Target']
            psm_df = psm_df[psm_df['pep_qvalue'] <= self.pep_fdr]
            if ignore_mod:
                peptide_df = psm_df.groupby('sequence', as_index=False).agg({
                    'peptide': lambda x: ','.join(set(x)),
                    'prev_AA': lambda x: ','.join(set(x)),
                    'next_AA': lambda x: ','.join(set(x)),
                    'seq_len': 'first',
                    'charge': lambda x: ','.join(sorted(set(x))),
                    'score': 'max',
                    'protein': lambda x: ','.join(set(p for ps in x for p in ps.split(';') \
                                                      if len(p.strip()) > 0 and not p.startswith(self.decoy_prefix)))
                })
                pep_stat = psm_df['sequence'].value_counts().reset_index().rename(columns={'count': file + '_spec_count'})
                peptide_df = peptide_df.merge(pep_stat, how='left', on='sequence')
            else:
                peptide_df = psm_df.groupby('peptide', as_index=False).agg({
                    'sequence': 'first',
                    'prev_AA': lambda x: ','.join(set(x)),
                    'next_AA': lambda x: ','.join(set(x)),
                    'seq_len': 'first',
                    'charge': lambda x: ','.join(sorted(set(x))),
                    'score': 'max',
                    'protein': lambda x: ','.join(set(p for ps in x for p in ps.split(';') \
                                                      if len(p.strip()) > 0 and not p.startswith(self.decoy_prefix)))
                })
                pep_stat = psm_df['peptide'].value_counts().reset_index().rename(columns={'count': file + '_spec_count'})
                peptide_df = peptide_df.merge(pep_stat, how='left', on='peptide')
            pep_result_map[file] = peptide_df
        return pep_result_map

    def combine_pep(self, pep_result_map, ignore_mod=False):
        common_columns = ['peptide', 'sequence', 'prev_AA', 'next_AA', 'seq_len', 'charge', 'score', 'protein']
        combined_pep_df = pd.DataFrame(columns=common_columns)

        def split_combine_unique(col1, col2):
            col1_split = col1.str.split(',').apply(lambda x: x if isinstance(x, list) else [])
            col2_split = col2.str.split(',').apply(lambda x: x if isinstance(x, list) else [])
            return (col1_split + col2_split).apply(lambda x: ','.join(sorted(set(x))))

        for pep_df in pep_result_map.values():
            if ignore_mod:
                combined_pep_df = combined_pep_df.merge(pep_df, on='sequence', how='outer', suffixes=('_df1', '_df2'))
                combined_pep_df['peptide'] = split_combine_unique(combined_pep_df['peptide_df1'], combined_pep_df['peptide_df2'])
            else:
                combined_pep_df = combined_pep_df.merge(pep_df, on='peptide', how='outer', suffixes=('_df1', '_df2'))
                combined_pep_df['sequence'] = split_combine_unique(combined_pep_df['sequence_df1'], combined_pep_df['sequence_df2'])

            combined_pep_df['prev_AA'] = split_combine_unique(combined_pep_df['prev_AA_df1'], combined_pep_df['prev_AA_df2'])
            combined_pep_df['next_AA'] = split_combine_unique(combined_pep_df['next_AA_df1'], combined_pep_df['next_AA_df2'])
            combined_pep_df['seq_len'] = split_combine_unique(combined_pep_df['seq_len_df1'], combined_pep_df['seq_len_df2'])
            combined_pep_df['charge'] = split_combine_unique(combined_pep_df['charge_df1'], combined_pep_df['charge_df2'])
            combined_pep_df['score'] =  combined_pep_df[['score_df1', 'score_df2']].max(axis=1)
            combined_pep_df['protein'] = split_combine_unique(combined_pep_df['protein_df1'], combined_pep_df['protein_df2'])
            combined_pep_df.drop(columns=[col for col in combined_pep_df.columns if col.endswith('_df1') or col.endswith('_df2')], inplace=True)

        spec_count_cols = [col for col in combined_pep_df.columns if col.endswith('_spec_count')]
        combined_pep_df[spec_count_cols] = combined_pep_df[spec_count_cols].fillna(0).astype(int)
        combined_pep_df = combined_pep_df[common_columns + spec_count_cols]
        return combined_pep_df

    def combine_protein(self):
        return

    def remove_contaminants(self, custom_fasta_path):
        return

if __name__ == '__main__':
    reporter = ResultCombiner('/mnt/d/workspace/mhc-booster/experiment/JY_Fractionation_Replicate_1/mhcbooster')
    reporter.run()
