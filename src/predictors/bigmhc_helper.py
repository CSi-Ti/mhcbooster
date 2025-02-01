
import tempfile
import subprocess
import time
import requests
import zipfile

import mhcgnomes
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List
from tqdm import tqdm
from mhcnames import normalize_allele_name
from src.utils.constants import EPSILON
from src.predictors.base_predictor_helper import BasePredictorHelper


class BigMhcHelper(BasePredictorHelper):
    def __init__(self,
                 peptides: list[str],
                 alleles: list[str],
                 report_directory: str):
        super().__init__('BigMHC', report_directory)

        self.bigmhc_url = 'https://github.com/KarchinLab/bigmhc/archive/refs/heads/master.zip'
        self.bigmhc_root = Path(__file__).parent.parent.parent / 'third_party' / 'bigmhc-master'
        if not self.bigmhc_root.exists():
            self._fetch_bigmhc_by_url()

        if alleles is None or len(alleles) == 0:
            raise RuntimeError('Alleles are needed for BigMHC predictions.')
        if np.min(np.vectorize(len)(peptides)) < 8:
            raise RuntimeError('BigMHC cannot make predictions on peptides shorter than 8 mer.')

        self.peptides = peptides
        self.alleles = self._format_class_I_alleles(alleles)
        self.bigmhc_exe_path = Path(__file__).parent.parent.parent / 'third_party' / 'bigmhc-master' / 'src' / 'predict.py'

    def _fetch_bigmhc_by_url(self):
        # Fetch BigMHC from GitHub
        third_party_root = self.bigmhc_root.parent
        if not third_party_root.exists():
            third_party_root.mkdir()

        print('Start to download BigMHC from GitHub. It will take several minutes...')
        response = requests.get(self.bigmhc_url, stream=True)
        if response.status_code == 200:
            destination_path = third_party_root / 'bigmhc-master.zip'
            with open(destination_path, 'wb') as file:
                # Initialize tqdm progress bar with total size and chunk size
                with tqdm(unit='B', unit_scale=True) as pbar:
                    # Download the file in chunks
                    for chunk in response.iter_content(chunk_size=1024):
                        # Write the chunk to the file
                        file.write(chunk)
                        # Update the progress bar with the size of the chunk
                        pbar.update(len(chunk))

                print("\nDownload completed!")
            zip_file = zipfile.ZipFile(str(destination_path))
            zip_file.extractall(third_party_root)
            print(f'Successfully downloaded BigMHC from GitHub and extracted to {third_party_root}')
        else:
            print(f'Failed to download BigMHC from GitHub. Status code: {response.status_code}')

    def _format_class_I_alleles(self, alleles: List[str]):
        std_alleles = []
        for allele in set(alleles):
            try:
                std_alleles.append(normalize_allele_name(allele))
            except ValueError:
                try:
                    std_alleles.append(mhcgnomes.parse(allele).to_string())
                except ValueError:
                    print(f'ERROR: Allele {allele} not supported.')
        return std_alleles

    def predict_df(self):
        print('Running BigMHC...')
        with tempfile.NamedTemporaryFile('w', delete=False) as pep_file:
            pep_file.write('mhc,pep\n')
            for pep in self.peptides:
                for allele in self.alleles:
                    pep_file.write(f'{allele},{pep}\n')
            pep_file_path = pep_file.name
        with tempfile.NamedTemporaryFile('w') as results:
            results_file_path = results.name

        command = f'python {self.bigmhc_exe_path} -i {pep_file_path} -o {results_file_path} -m el -d cpu'
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        for line in iter(process.stdout.readline, ''):
            print(line.strip())
        for line in iter(process.stderr.readline, ''):
            print(line.strip())
        process.stdout.close()
        process.stderr.close()
        process.wait()

        self.pred_df = pd.read_csv(results_file_path, index_col=False)

        return self.pred_df

    def score_df(self) -> pd.DataFrame:
        """
        Add features from mhcflurry predictions to feature_matrix. Affinity predictions added as log values.
        All non-log value clipped to a minimum of 1e-7.
        """

        predictions = pd.DataFrame()

        alleles = list(self.pred_df.loc[:, 'mhc'].unique())
        for allele in alleles:
            df = self.pred_df.loc[self.pred_df['mhc'] == allele, :]
            assert list(df['pep']) == list(self.peptides)
            predictions[f'{allele}_BigMHC_ELScore'] = df['BigMHC_EL'].clip(lower=EPSILON).to_numpy()
            predictions[f'{allele}_logBigMHC_ELScore'] = np.log(df['BigMHC_EL'].clip(lower=EPSILON)).to_numpy()
        return predictions
