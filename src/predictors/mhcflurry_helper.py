
import tempfile
import subprocess

import numpy as np
import pandas as pd
from typing import List
from mhcnames import normalize_allele_name

from mhcvalidator.constants import EPSILON
from mhcvalidator.predictors.base_predictor_helper import BasePredictorHelper


class MhcFlurryHelper(BasePredictorHelper):
    def __init__(self,
                 peptides: list[str],
                 alleles: list[str],
                 report_directory: str):
        super().__init__('MHCflurry', report_directory)

        if alleles is None or len(alleles) == 0:
            raise RuntimeError('Alleles are needed for MhcFlurry predictions.')
        if np.max(np.vectorize(len)(peptides)) > 16:
            raise RuntimeError('MhcFlurry cannot make predictions on peptides over length 16.')

        self.peptides = peptides
        self.alleles = self._format_class_I_alleles(alleles)

    def _format_class_I_alleles(self, alleles: List[str]):
        std_alleles = []
        for allele in set(alleles):
            try:
                std_alleles.append(normalize_allele_name(allele))
            except ValueError:
                print(f'ERROR: Allele {allele} not supported.')
        return [a.replace('*', '').replace(':', '') for a in std_alleles]

    def predict_df(self):
        # we will run MhcFlurry in a separate process so the Tensorflow space doesn't get messed up. I don't know why it
        # happens, but it does, and then either MhcFlurry or TensorFlow Decision Forests stops working.
        # I think perhaps because MhcFlurry uses some legacy code from TFv1 (I think), though this is only
        # a suspicion.
        print('Running MhcFlurry')
        with tempfile.NamedTemporaryFile('w', delete=False) as pep_file:
            pep_file.write('allele,peptide\n')
            for pep in self.peptides:
                for allele in self.alleles:
                    pep_file.write(f'{allele},{pep}\n')
            pep_file_path = pep_file.name
        with tempfile.NamedTemporaryFile('w') as results:
            results_file = results.name

        command = f'mhcflurry-predict --out {results_file} {pep_file_path}'.split()
        p = subprocess.Popen(command)
        _ = p.communicate()

        self.pred_df = pd.read_csv(results_file, index_col=False)

        return self.pred_df

    def score_df(self) -> pd.DataFrame:
        """
        Add features from mhcflurry predictions to feature_matrix. Affinity predictions added as log values.
        All non-log value clipped to a minimum of 1e-7.
        """

        predictions = pd.DataFrame()

        alleles = list(self.pred_df.loc[:, 'allele'].unique())
        for allele in alleles:
            df = self.pred_df.loc[self.pred_df['allele'] == allele, :]
            assert list(df['peptide']) == list(self.peptides)
            predictions[f'{allele}_MhcFlurry_PresentationScore'] = df['mhcflurry_presentation_score'].clip(lower=EPSILON).to_numpy()
            predictions[f'{allele}_logMhcFlurry_PresentationScore'] = np.log(predictions[f'{allele}_MhcFlurry_PresentationScore'] + 0.01)
            predictions[f'{allele}_logMhcFlurry_Affinity'] = np.log(df['mhcflurry_affinity'].clip(lower=EPSILON)).to_numpy()
        df = self.pred_df.loc[self.pred_df['allele'] == alleles[0], :]
        # we only need one processing score column
        predictions[f'MhcFlurry_ProcessingScore'] = df['mhcflurry_processing_score'].clip(lower=EPSILON).to_numpy()
        predictions[f'logMhcFlurry_ProcessingScore'] = np.log(predictions[f'MhcFlurry_ProcessingScore'] + 0.01)
        return predictions

    def format_pred_result_for_saving(self) -> pd.DataFrame:
        predictions = pd.DataFrame()
        predictions['Peptide'] = np.array(self.peptides)

        alleles = list(self.pred_df.loc[:, 'allele'].unique())
        for allele in alleles:
            df = self.pred_df.loc[self.pred_df['allele'] == allele, :]
            assert list(df['peptide']) == list(self.peptides)
            for pred in df.columns:
                if pred in ['allele', 'peptide', 'mhcflurry_processing_score']:
                    continue
                predictions[f'{allele}_MhcFlurry_{pred}'] = df[pred].clip(lower=EPSILON).to_numpy()
            predictions[f'{allele}_logMhcFlurry_Affinity'] = np.log(df['mhcflurry_affinity'].clip(lower=EPSILON)).to_numpy()
        df = self.pred_df.loc[self.pred_df['allele'] == alleles[0], :]
        predictions[f'MhcFlurry_ProcessingScore'] = df['mhcflurry_processing_score'].clip(lower=EPSILON).to_numpy()
        return predictions