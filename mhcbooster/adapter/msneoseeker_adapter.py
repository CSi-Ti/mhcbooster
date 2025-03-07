import argparse

import numpy as np
import pandas as pd
from pathlib import Path
from mhcnames import normalize_allele_name


description = 'MHCBooster result converter for MSNeoSeeker'
parser = argparse.ArgumentParser(description=description)
parser.add_argument('--result_dir',
                     required=True,
                     type=str,
                     help='Result directory of MHCBooster to convert.')
parser.add_argument('--output_dir',
                    required=False,
                    type=str,
                    help='Output directory for converted files.')

def convert_protein_id(df):
    assert 'protein_id' in df.columns
    id_col_idx = df.columns.get_loc('protein_id')
    protein_ids = df['protein_id'].apply(lambda x: str(x).split('-')[0] if '-' in str(x) else str(x))
    transcript_ids = df['protein_id'].apply(lambda x: str(x).split('-')[1] if '-' in str(x) else '')
    df['protein_id'] = protein_ids
    if 'transcript_id' not in df.columns:
        df.insert(id_col_idx + 1, 'transcript_id', transcript_ids)

def prepare_app_df(sequence, result_dir):
    app_df = pd.DataFrame()
    app_df['sequence'] = sequence
    app_score_paths = Path(result_dir).glob('app_prediction.*.tsv')
    for app_score_path in app_score_paths:
        predictor = app_score_path.stem.split('.')[1]
        psm_app_df = pd.read_csv(app_score_path, sep='\t')
        if predictor == 'netmhcpan':
            seq_app_df = psm_app_df.loc[
                psm_app_df.groupby('Peptide')['EL_Rank'].idxmin(), ['Peptide', 'Allele', 'EL_Rank']]
            seq_app_df['best_predictor'] = predictor
            seq_app_df['netmhcpan_binder'] = seq_app_df['EL_Rank'].apply(
                lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            seq_app_df['netmhcpan_rank'] = seq_app_df['EL_Rank']
            seq_app_df['netmhcpan_best_allele'] = seq_app_df['Allele']
        if predictor == 'mhcflurry':
            seq_app_df = psm_app_df.loc[
                psm_app_df.groupby('peptide')['mhcflurry_affinity_percentile'].idxmin(), ['peptide', 'allele',
                                                                                          'mhcflurry_affinity_percentile']]
            seq_app_df['best_predictor'] = predictor
            seq_app_df['mhcflurry_binder'] = seq_app_df['mhcflurry_affinity_percentile'].apply(
                lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            seq_app_df['mhcflurry_rank'] = seq_app_df['mhcflurry_affinity_percentile']
            seq_app_df['mhcflurry_best_allele'] = seq_app_df['allele']
        if predictor == 'bigmhc':
            seq_app_df = psm_app_df.loc[psm_app_df.groupby('pep')['BigMHC_EL'].idxmin(), ['pep', 'mhc', 'BigMHC_EL']]
            seq_app_df['best_predictor'] = predictor
            seq_app_df['BigMHC_EL'] *= 100
            seq_app_df['bigmhc_binder'] = seq_app_df['BigMHC_EL'].apply(
                lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            seq_app_df['bigmhc_rank'] = seq_app_df['BigMHC_EL']
            seq_app_df['bigmhc_best_allele'] = seq_app_df['mhc']
        if predictor == 'netmhciipan':
            seq_app_df = psm_app_df.loc[
                psm_app_df.groupby('Peptide')['EL_Rank'].idxmin(), ['Peptide', 'Allele', 'EL_Rank']]
            seq_app_df['best_predictor'] = predictor
            seq_app_df['netmhciipan_binder'] = seq_app_df['EL_Rank'].apply(
                lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            seq_app_df['netmhciipan_rank'] = seq_app_df['EL_Rank']
            seq_app_df['netmhciipan_best_allele'] = seq_app_df['Allele']
        if predictor == 'mixmhc2pred':
            psm_app_df['BestAllele'] = psm_app_df['BestAllele'].fillna('')
            psm_app_df['%Rank_best'] = psm_app_df['%Rank_best'].fillna(100)
            seq_app_df = psm_app_df.loc[
                psm_app_df.groupby('Peptide')['%Rank_best'].idxmin(), ['Peptide', 'BestAllele', '%Rank_best']]
            seq_app_df['best_predictor'] = predictor
            seq_app_df['mixmhc2pred_binder'] = seq_app_df['%Rank_best'].apply(
                lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            seq_app_df['mixmhc2pred_rank'] = seq_app_df['%Rank_best']
            seq_app_df['mixmhc2pred_best_allele'] = seq_app_df['BestAllele']
        seq_app_df.columns = ['sequence', 'best_allele', 'min_rank', 'best_predictor'] + [col for col in seq_app_df.columns[4:]]
        app_df = pd.merge(app_df, seq_app_df, on='sequence', how='left', suffixes=('_left', '_right'))
        if 'min_rank_left' in app_df.columns:
            app_df['min_rank'] = app_df[['min_rank_left', 'min_rank_right']].min(axis=1)
            app_df['best_allele'] = np.where(app_df['min_rank'] == app_df['min_rank_left'], app_df['best_allele_left'],
                                             app_df['best_allele_right'])
            app_df['best_predictor'] = np.where(app_df['min_rank'] == app_df['min_rank_left'], app_df['best_predictor_left'],
                                             app_df['best_predictor_right'])
            app_df.drop(columns=['min_rank_left', 'min_rank_right', 'best_allele_left', 'best_allele_right',
                                 'best_predictor_left', 'best_predictor_right'],
                        inplace=True)
        for col in app_df.columns:
            if 'allele' in col:
                app_df[col] = app_df[col].apply(lambda a: normalize_allele_name(a.replace('__', '-')))
    return app_df

def convert_app_info(df, result_dir):
    app_df = prepare_app_df(df['sequence'], result_dir=result_dir)
    binder_col_idx = [i for i in range(len(df.columns)) if 'binder' in df.columns[i]]
    if len(binder_col_idx) > 0:
        col_to_remove = [df.columns[i] for i in range(binder_col_idx[0], binder_col_idx[-1] + 1)]
        df.drop(columns=col_to_remove, inplace=True)
        for i in range(1, len(app_df.columns)):
            idx_to_insert = binder_col_idx[0] + i - 1
            col_name = app_df.columns[i]
            df.insert(idx_to_insert, col_name, app_df[col_name])


def run(result_dir, output_dir):

    psm_paths = Path(result_dir).rglob('psm.tsv')
    for psm_path in psm_paths:
        output_path = Path(output_dir)/psm_path.parent.name
        output_path.mkdir(parents=True, exist_ok=True)

        # psm
        psm_df = pd.read_csv(psm_path, sep='\t')
        convert_protein_id(psm_df)
        convert_app_info(psm_df, psm_path.parent)
        psm_df.to_csv(output_path/'psm.msneoseeker.tsv', sep='\t', index=False)

        # peptide
        peptide_path = psm_path.parent/'peptide.tsv'
        peptide_df = pd.read_csv(peptide_path, sep='\t')
        convert_protein_id(peptide_df)
        convert_app_info(peptide_df, peptide_path.parent)
        peptide_df.to_csv(output_path/'peptide.msneoseeker.tsv', sep='\t', index=False)

        # sequence
        sequence_path = psm_path.parent/'sequence.tsv'
        sequence_df = pd.read_csv(sequence_path, sep='\t')
        convert_protein_id(sequence_df)
        convert_app_info(sequence_df, peptide_path.parent)
        sequence_df.to_csv(output_path/'sequence.msneoseeker.tsv', sep='\t', index=False)


if __name__ == "__main__":
    args = parser.parse_args()
    assert args.result_dir != args.output_dir, 'Output directory cannot be the same as result directory.'

    if args.output_dir is None:
        args.output_dir = Path(args.result_dir) / 'msneoseeker_format'

    run(args.result_dir, args.output_dir)
