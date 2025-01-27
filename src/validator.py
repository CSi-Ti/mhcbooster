import os
import re
import time

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import pandas as pd
import numpy as np
from numpy.random import RandomState
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import backend as K
from typing import Union, List
from os import PathLike
from pathlib import Path
from mhcvalidator.data_loaders import load_file, load_pout_data
from mhcvalidator.tools.mzml_parser import get_rt_ccs_ms2_from_mzml, get_rt_ccs_ms2_from_msfragger_mzml
from mhcvalidator.features import prepare_features
from mhcvalidator.predictors.netmhcpan_helper import NetMHCpanHelper
from mhcvalidator.predictors.mhcflurry_helper import MhcFlurryHelper
from mhcvalidator.predictors.bigmhc_helper import BigMhcHelper
from mhcvalidator.predictors.mixmhc2pred_helper import MixMhc2PredHelper
from mhcvalidator.predictors.peptdeep_helper import PeptDeepHelper
from mhcvalidator.predictors.autort_helper import AutortHelper
from mhcvalidator.predictors.deeplc_helper import DeepLCHelper
from mhcvalidator.predictors.im2deep_helper import IM2DeepHelper
from mhcvalidator.predictors.koina_helper import KoinaHelper, KOINA_PREDICTORS
from mhcvalidator.fdr import calculate_qs, calculate_peptide_level_qs, calculate_roc
import matplotlib.pyplot as plt
from mhcflurry.encodable_sequences import EncodableSequences
from mhcvalidator.models import get_model_without_peptide_encoding, get_model_with_peptide_encoding
from mhcvalidator.peptides import remove_previous_and_next_aa, clean_peptide_sequences, remove_charge, \
    remove_modifications
from mhcvalidator.nd_standard_scalar import NDStandardScaler
from copy import deepcopy
from datetime import datetime
import tempfile
from collections import Counter
import tensorflow.python.util.deprecation as deprecation
from matplotlib.gridspec import GridSpec
import matplotlib.backends.backend_pdf as plt_pdf
from inspect import signature
from matplotlib.cm import get_cmap
from mhcvalidator.datasets import k_fold_split, k_fold_split_s, k_fold_split_so

deprecation._PRINT_DEPRECATION_WARNINGS = False

# This can be uncommented to prevent the GPU from getting used.
# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

DEFAULT_TEMP_MODEL_DIR = str(Path(tempfile.gettempdir()) / 'validator_models')


class MhcValidator:
    def __init__(self,
                 random_seed: int = 0,
                 model_dir: Union[str, PathLike] = DEFAULT_TEMP_MODEL_DIR,
                 max_threads: int = -1):
        self.filename: Union[str, None] = None
        self.filepath: Union[Path, None] = None
        self.model: keras.Model = None
        self.raw_data: Union[pd.DataFrame, None] = None
        self.feature_matrix: Union[pd.DataFrame, None] = None
        self.labels: Union[List[int], None] = None
        self.peptides: Union[List[str], None] = None
        self.peptides_with_mods: Union[List[str], None] = None
        self.charges: Union[List[int], None] = None
        self.high_prob_indices: Union[np.ndarray, None] = None
        self.exp_rts: Union[np.ndarray, None] = None
        self.exp_ims: Union[np.ndarray, None] = None
        self.exp_spectra: Union[pd.DataFrame, None] = None
        self.encoded_peptides = None
        self.loaded_filetype: Union[str, None] = None
        self.random_seed: int = random_seed
        self.predictions: np.array = None
        self.qs: np.array = None
        self.roc = None
        self.percolator_qs = None
        self.obs_rts: np.array = None
        self.pred_rts: np.array = None
        self.mhc_class: Union[str, None] = None
        self.alleles: List[str] = None
        self.min_len: int = 5
        self.max_len: int = 100
        self.model_dir = Path(model_dir)
        if not os.path.exists(self.model_dir):
            os.mkdir(self.model_dir)
        self.fine_tune: bool = False
        self.koina_predictors: List[str] = []
        self.mhcflurry_predictions: pd.DataFrame = None
        self.netmhcpan_predictions: pd.DataFrame = None
        self.peptdeep_predictions: pd.DataFrame = None
        self.annotated_data: pd.DataFrame = None
        if max_threads < 1:
            self.max_threads: int = os.cpu_count()
        else:
            self.max_threads: int = max_threads

    def set_mhc_params(self,
                       alleles: Union[str, List[str]] = None,
                       mhc_class: str = 'I',
                       max_pep_len: int = None,
                       min_pep_len: int = None) -> None:
        """
        Set the MHC-specific parameters.

        :param alleles: The alleles to be used by MhcFlurry or NetMHCpan.
        :param mhc_class: The MHC class of the peptides. Must be one of {'I', 'II'}
        :param min_pep_len: Maximum length of peptides allowed. Will default to 16 for class I and 30 for class II. Note
        that MhcFlurry does not accept peptide lengths greater than 16. There is no length restriction for NetMHCpan.
        :param max_pep_len: Minimum length of peptides allowed. Will default to 8 for class I and 9 for class II. Note
        that NetMHC(II)pan does not accept peptide lengths less than 8 for class I or 9 for class I. NetMHCpan predictions
        take much longer for longer peptides.
        :return: None
        """

        if isinstance(alleles, str):
            alleles = [alleles]
        self.alleles = list(set(alleles))

        if min_pep_len is not None:
            self.min_len = min_pep_len
        if max_pep_len is not None:
            self.max_len = max_pep_len

        if mhc_class is not None:
            if mhc_class not in ['I', 'II']:
                raise ValueError("mhc_class must be one of {'I', 'II'}")
            self.mhc_class = mhc_class

        if mhc_class == 'I' and self.min_len < 8:
            self.min_len = 8
            print('The minimum peptide length is reset to 8, due to limitations in MHC-I predictors.')
        if mhc_class == 'I' and self.max_len > 15:
            self.max_len = 15
            print('The maximum peptide length is reset to 15, due to limitations in MHC-I predictors.')
        if mhc_class == 'II' and self.min_len < 9:
            self.min_len = 9
            print('The minimum peptide length is reset to 9, due to limitations in MHC-II predictors.')
        if mhc_class == 'II' and self.max_len > 30:
            self.max_len = 30
            print('The maximum peptide length is reset to 30, due to limitations in MHC-II predictors.')


    def _check_peptide_lengths(self):
        max_len = self.max_len
        longest_peptide = np.max(np.vectorize(len)(self.peptides))
        if max_len > longest_peptide:
            print(f'Longest peptide ({longest_peptide} mer) is shorter than set maximum length ({max_len} mer). '
                  f'Changing max_len to {longest_peptide}.')
            self.max_len = longest_peptide

    def load_data(self,
                  filepath: Union[str, PathLike],
                  filetype='auto',
                  decoy_tag='rev_',
                  peptide_column: str = None,
                  protein_column: str = None,
                  tag_is_prefix: bool = True,
                  file_delimiter: str = '\t',
                  use_features: Union[List[str], None] = None):
        """
        Load the results of an upstream search or validation tool. PIN, pepXML, mzID, X! Tandem, Spectromine and
        generic tabular formats are accepted. To load POUT files, use the separate 'load_pout_data' function. You can
        load both PIN and POUT files from a single experiment using the separate 'load_percolator_data' function.
        Generic tabular files must contain a column titled 'Peptide' or 'peptide' which contains the peptide sequences.

        :param filepath: The path to the file you want to load. Can be absolute or relative.
        :param filetype: The type of file. Must be one of {'auto', 'pin', 'pepxml', 'tabular', 'mzid', 'tandem',
            'spectromine'}. If you choose 'auto', the file type will be inferred from the file extension. Be
            cautious when loading pepXML and X! Tandem files, as the extensions are similar. It is best to be explicit
            in these cases.
        :param decoy_tag: The decoy tag used in the upstream FASTA to differentiate decoy proteins from target proteins.
        :param protein_column: The header of the column containing the protein IDs. Required for tabular data of an
            unspecified format.
        :param tag_is_prefix: Whether or not the decoy tag is a prefix. If False, it is assumed the tag is a suffix.
        :param file_delimiter: The delimiter used if the file is tabular.
        :param use_features: A list of column headers to be used as training features. Not required  If your tabular data
            contains a column indicating the target/decoy label of each PSM, DO NOT INCLUDE THIS COLUMN! The label will
            be determined from the protein IDs.
        :return: None
        """

        if filetype == 'auto':
            if str(filepath).lower().endswith('pin') or str(filepath).lower().endswith('mhcv'):
                filetype = 'pin'
            elif str(filepath).lower().endswith('pepxml'):
                filetype = 'tandem'
            elif str(filepath).lower().endswith('pep.xml'):
                filetype = 'pepxml'
            else:
                raise ValueError('File type could not be inferred from filename. You must explicitly specify the '
                                 'filetype.')
        else:
            if filetype not in ['auto', 'pin', 'pepxml', 'tabular', 'mhcv']:
                raise ValueError("filetype must be one of "
                                 "{'auto', 'pin', 'pepxml', 'tabular', 'mhcv'}")

        print(f'MHC class: {self.mhc_class if self.mhc_class else "not specified"}')
        print(f'Alleles: {self.alleles if self.alleles else "not specified"}')
        print(f'Minimum peptide length: {self.min_len}')
        print(f'Maximum peptide length: {self.max_len}')

        print('Loading PSM file...')
        self.raw_data = load_file(filename=filepath, filetype=filetype, decoy_tag=decoy_tag,
                                  protein_column=protein_column, file_sep=file_delimiter,
                                  tag_is_prefix=tag_is_prefix, min_len=self.min_len, max_len=self.max_len)
        self.labels = self.raw_data['Label'].to_numpy()

        # Peptide
        if peptide_column is not None:
            self.peptides_with_mods = list(self.raw_data[peptide_column])
        elif filetype == 'pin':
            self.peptides_with_mods = list(self.raw_data['Peptide'])
        #elif filetype == 'mzid':
        #    self.peptides = list(self.raw_data['PeptideSequence'])
        #elif filetype == 'spectromine':
        #    self.peptides = list(self.raw_data['PEP.StrippedSequence'])
        else:
            if 'peptide' in self.raw_data.columns:
                self.peptides_with_mods = list(self.raw_data['peptide'])
            elif 'Peptide' in self.raw_data.columns:
                self.peptides_with_mods = list(self.raw_data['Peptide'])
            else:
                raise IndexError('Peptide field could not be automatically found. Please indicate the column '
                                 'containing the peptide sequences')
        self.peptides_with_mods = remove_charge(remove_previous_and_next_aa(self.peptides_with_mods))
        self.peptides = np.array(remove_modifications(self.peptides_with_mods))
        self._check_peptide_lengths()

        # Charge
        self.charges = np.zeros(len(self.peptides_with_mods), dtype=int)
        for col in self.raw_data.columns:
            if 'charge' in col.lower():
                self.charges[self.raw_data[col] == '1'] =  int(re.findall(r'\d+', col)[0])

        # High prob indices
        self.high_prob_indices = None
        qs_threshold = 0.0001
        max_qs_threshold = 0.01
        min_points = 100

        if 'lnExpect' in self.raw_data.columns:
            qs = calculate_qs(self.raw_data['lnExpect'].astype(float), self.labels, higher_better=False)
        elif 'log10_evalue' in self.raw_data.columns:
            qs = calculate_qs(self.raw_data['log10_evalue'].astype(float), self.labels, higher_better=False)
        else:
            qs = None
            print('lnExpect or log10_evalue score cannot be found from input files. Processing without calibration!')

        if qs is not None:
            high_prob_indices = qs < qs_threshold
            if np.sum(high_prob_indices) >= min_points:
                self.high_prob_indices = high_prob_indices
            else:
                tmp_qs_threshold = np.sort(qs)[min_points]
                if tmp_qs_threshold < max_qs_threshold:
                    self.high_prob_indices = qs < tmp_qs_threshold
                    print(f'Not enough PSMs for calibration. Relaxed the high confidence q-value threshold to {tmp_qs_threshold}')
                else:
                    print('Not enough PSMs for calibration. Processing without calibration!')

        print(f'Loaded {len(self.peptides)} PSMs, including {np.sum(self.high_prob_indices)} high confidence PSMs.')

        self.loaded_filetype = filetype
        self.filename = Path(filepath).name
        self.filepath = Path(filepath).expanduser().resolve()

        print('Preparaing training features')
        self.feature_matrix = prepare_features(self.raw_data, filetype=self.loaded_filetype, use_features=use_features)

    def load_pout_data(self,
                       targets_pout: Union[str, PathLike],
                       decoys_pout: Union[str, PathLike],
                       use_features: Union[List[str], None] = None) -> None:
        """
        Load POUT files generated by Percolator. You must have created both target and decoy POUT files from Percolator.

        :param targets_pout: The path to the targets POUT file.
        :param decoys_pout: The path to the decoys POUT file.
        :param use_features: (Optional) A list of features (i.e. columns) to load.
        :return: None
        """

        print(f'MHC class: {self.mhc_class if self.mhc_class else "not specified"}')
        print(f'Alleles: {self.alleles if self.alleles else "not specified"}')
        print(f'Minimum peptide length: {self.min_len}')
        print(f'Maximum peptide length: {self.max_len}')

        print('Loading PSM file')
        self.raw_data = load_pout_data(targets_pout, decoys_pout, self.min_len, self.max_len)
        self.labels = self.raw_data['Label'].values
        self.peptides = list(self.raw_data['peptide'])
        self.peptides = np.array(clean_peptide_sequences(self.peptides))
        self._check_peptide_lengths()

        # self.raw_data.drop(columns=['Label'], inplace=True)
        self.loaded_filetype = 'tabular'
        self.filename = (Path(targets_pout).name, Path(decoys_pout).name)
        self.filepath = (Path(targets_pout).expanduser().resolve(), Path(decoys_pout).expanduser().resolve())

        self.feature_matrix = prepare_features(self.raw_data,
                                               filetype=self.loaded_filetype,
                                               use_features=use_features)

    def load_percolator_data(self,
                             pin_file: Union[str, PathLike],
                             target_pout_file: Union[str, PathLike],
                             decoy_pout_file: Union[str, PathLike],
                             use_features: Union[List[str], None] = None,
                             decoy_tag='rev_',
                             tag_is_prefix: bool = True
                             ) -> None:
        """
        Load PIN and POUT files from a single experiment.You must have created both target and decoy POUT files from
        Percolator.

        :param pin_file: Path to the PIN file.
        :param target_pout_file: The path to the targets POUT file.
        :param decoy_pout_file: The path to the decoys POUT file.
        :param use_features: (Optional) A list of features (i.e. columns) to load.
        :param decoy_tag: The decoy tag used to indicate decoys in the upstream FASTA file.
        :param tag_is_prefix: Whether or not the decoy tag is a prefix. If False, it is assumed the tag is a suffix.
        :return: None
        """
        print(f'MHC class: {self.mhc_class if self.mhc_class else "not specified"}')
        print(f'Alleles: {self.alleles if self.alleles else "not specified"}')
        print(f'Minimum peptide length: {self.min_len}')
        print(f'Maximum peptide length: {self.max_len}')
        print('Loading PSM file')

        pout_data = load_pout_data(target_pout_file, decoy_pout_file, self.min_len, self.max_len)

        pin_data = load_file(filename=pin_file, filetype='pin', decoy_tag=decoy_tag,
                             protein_column='Proteins', file_sep='\t',
                             tag_is_prefix=tag_is_prefix, min_len=self.min_len, max_len=self.max_len)
        pout_data.drop(columns=['peptide', 'proteinIds'], inplace=True)
        pin_data.drop(columns=['Label'], inplace=True)

        self.raw_data = pin_data.join(pout_data.set_index('PSMId'), on='SpecId')
        self.percolator_qs = self.raw_data['q-value'].to_numpy(np.float32)
        self.raw_data.drop(columns=['q-value'], inplace=True)

        self.labels = self.raw_data['Label'].to_numpy(np.float32)
        self.peptides = list(self.raw_data['Peptide'])
        self.peptides = np.array(clean_peptide_sequences(self.peptides))
        self._check_peptide_lengths()

        self.loaded_filetype = 'PIN_POUT'
        self.filename = (Path(pin_file).name,
                         Path(target_pout_file).name,
                         Path(target_pout_file).name)
        self.filepath = (Path(pin_file).expanduser().resolve(),
                         Path(target_pout_file).expanduser().resolve(),
                         Path(target_pout_file).expanduser().resolve())
        self.feature_matrix = prepare_features(self.raw_data,
                                               filetype='pin',  # PIN file processing works for this
                                               use_features=use_features)

    def encode_peptide_sequences(self):
        """
        Use a BLOSUM62 substitution matrix to numerically encode each peptide sequence. Uses the EncodableSequences
        class from MhcFlurry. Encoded peptides are saved in self.encoded_peptides.
        :return:
        """

        encoder = EncodableSequences(list(self.peptides))
        padding = 'pad_middle' if self.mhc_class == 'I' else 'left_pad_right_pad'
        encoded_peps = encoder.variable_length_to_fixed_length_vector_encoding('BLOSUM62',
                                                                               max_length=self.max_len,
                                                                               alignment_method=padding)
        self.encoded_peptides = deepcopy(encoded_peps)

    def add_mhcflurry_predictions(self):
        """
        Run MhcFlurry and add presentation predictions to the training feature matrix.

        :return: None
        """
        if self.mhc_class == 'II':
            raise RuntimeError('MhcFlurry is only compatible with MHC class I')

        mhcflurry_helper = MhcFlurryHelper(self.peptides, self.alleles, self.report_directory)
        mhcflurry_helper.predict_df()
        predictions = mhcflurry_helper.score_df()
        self.feature_matrix = self.feature_matrix.join(predictions)
        mhcflurry_helper.draw_prediction_distributions(predictions, self.labels)
        self.mhcflurry_predictions = mhcflurry_helper.format_pred_result_for_saving()


    def add_netmhcpan_predictions(self):
        """
        Run NetMHCpan and add presentation predictions to the training feature matrix.

        :return: None
        """
        if self.alleles is None or self.mhc_class is None:
            raise RuntimeError('You must first set the MHC parameters using mhcvalidator.set_mhc_params')

        print(f'Running NetMHC{"II" if self.mhc_class == "II" else ""}pan')
        netmhcpan_helper = NetMHCpanHelper(peptides=self.peptides, alleles=self.alleles,
                                           mhc_class=self.mhc_class, n_threads=self.max_threads,
                                           report_directory=self.report_directory)
        netmhcpan_helper.predict_df()
        predictions = netmhcpan_helper.score_df()
        self.feature_matrix = self.feature_matrix.join(predictions)
        netmhcpan_helper.draw_prediction_distributions(predictions, self.labels)
        self.netmhcpan_predictions = netmhcpan_helper.format_pred_result_for_saving()


    def add_bigmhc_predictions(self):
        if self.alleles is None or self.mhc_class is None:
            raise RuntimeError('You must first set the MHC parameters using mhcvalidator.set_mhc_params')

        bigmhc_helper = BigMhcHelper(peptides=self.peptides, alleles=self.alleles, report_directory=self.report_directory)
        bigmhc_helper.predict_df()
        predictions = bigmhc_helper.score_df()
        self.feature_matrix = self.feature_matrix.join(predictions)
        bigmhc_helper.draw_prediction_distributions(predictions, self.labels)


    def add_mixmhc2pred_predictions(self):
        if self.alleles is None:
            raise RuntimeError('You must first set the MHC parameters using mhcvalidator.set_mhc_params')

        mixmhc2pred_helper = MixMhc2PredHelper(peptides=self.peptides, alleles=self.alleles, report_directory=self.report_directory)
        mixmhc2pred_helper.predict_df()
        predictions = mixmhc2pred_helper.score_df()
        self.feature_matrix = self.feature_matrix.join(predictions)
        mixmhc2pred_helper.draw_prediction_distributions(predictions, self.labels)


    def add_autort_predictions(self):
        """
        Run AutoRT and add predicted RT scores to the training feature matrix.

        :return: None
        """
        autort_helper = AutortHelper(self.peptides_with_mods, self.exp_rts, self.high_prob_indices, fine_tune=self.fine_tune, verbose=False, report_directory=self.report_directory)
        autort_helper.predict_df()
        predictions = autort_helper.score_df()
        self.feature_matrix = self.feature_matrix.join(predictions)
        autort_helper.draw_prediction_distributions(predictions, self.labels)


    def add_deeplc_predictions(self):
        """
        Run DeepLC and add predicted RT scores to the training feature matrix.

        :return: None
        """
        deeplc_helper = DeepLCHelper(self.peptides, self.peptides_with_mods, self.exp_rts, self.high_prob_indices, fine_tune=self.fine_tune, verbose=False, report_directory=self.report_directory)
        deeplc_helper.predict_df()
        predictions = deeplc_helper.score_df()
        self.feature_matrix = self.feature_matrix.join(predictions)
        deeplc_helper.draw_prediction_distributions(predictions, self.labels)


    def add_im2deep_predictions(self):
        """
        Run IM2Deep and add predicted CCS scores to the training feature matrix.

        :return: None
        """
        if np.max(self.exp_ims) == 0:
            print('Cannot read ion mobility from experimental data. Skipping IM2Deep predictions...')
            return
        im2deep_helper = IM2DeepHelper(self.peptides, self.peptides_with_mods, self.charges, self.exp_ims,
                                       self.high_prob_indices, report_directory=self.report_directory,
                                       fine_tune=self.fine_tune, verbose=False)
        im2deep_helper.predict_df()
        predictions = im2deep_helper.score_df()
        self.feature_matrix = self.feature_matrix.join(predictions)
        im2deep_helper.draw_prediction_distributions(predictions, self.labels)


    def add_peptdeep_predictions(self):
        """
        Run AlphaPeptDeep and add presentation predictions to the training feature matrix.

        :return: None
        """
        peptdeep_helper = PeptDeepHelper(self.peptides, self.raw_data, self.report_directory)
        peptdeep_helper.predict_df()
        predictions = peptdeep_helper.score_df()
        self.feature_matrix = self.feature_matrix.join(predictions)
        peptdeep_helper.draw_prediction_distributions(predictions, self.labels)


    def add_koina_predictions(self):
        attempt = 0
        max_retries = 5  # Maximum number of retries

        while attempt < max_retries:
            try:
                koina_helper = KoinaHelper(self.peptides, self.peptides_with_mods, self.charges,
                                           predictor_names=self.koina_predictors,
                                           exp_rts=self.exp_rts, exp_ims=self.exp_ims, exp_spectra=self.exp_spectra,
                                           high_prob_indices=self.high_prob_indices,
                                           instrument_type='QE',
                                           fragmentation_type='HCD',
                                           koina_server_url=self.koina_server_url,
                                           report_directory=self.report_directory
                                           )
                koina_helper.predict_df()
                predictions = koina_helper.score_df()
                self.feature_matrix = self.feature_matrix.join(predictions)
                koina_helper.draw_prediction_distributions(predictions, self.labels)
                break
            except Exception as e:
                print(f"Warning: {e}")
                attempt += 1
                if attempt < max_retries:
                    print(f"Retrying in 1 minute... (Attempt {attempt}/{max_retries})")
                    time.sleep(60)
                else:
                    print("Max retries reached. Operation failed.")
                    raise


    def _set_seed(self, random_seed: int = None):
        if random_seed is None:
            random_seed = self.random_seed
        tf.random.set_seed(random_seed)
        np.random.seed(random_seed)


    def get_nn_model(self,
                     learning_rate: float = 0.001,
                     dropout: float = 0.5,
                     hidden_layers: int = 2,
                     width_ratio: float = 5.0,
                     loss_fn=tf.losses.BinaryCrossentropy()
                     ):
        """
        Return a compiled multilayer perceptron neural network with the indicated architecture.

        :param learning_rate: Learning rate used by the optimizer (adam).
        :param dropout: Dropout between each layer.
        :param hidden_layers: Number of hidden layers.
        :param width_ratio: Ratio of width of hidden layers to width of input layer.
        :param loss_fn: The loss function to use.
        :return: A compiled keras.Model
        """

        optimizer = keras.optimizers.Adam(learning_rate=learning_rate)
        model = get_model_without_peptide_encoding(self.feature_matrix.shape[1],
                                                   dropout=dropout,
                                                   hidden_layers=hidden_layers,
                                                   max_pep_length=self.max_len,
                                                   width_ratio=width_ratio)
        model.compile(loss=loss_fn, optimizer=optimizer)

        return model

    def get_nn_model_with_sequence_encoding(self,
                                            learning_rate: float = 0.001,
                                            dropout: float = 0.5,
                                            hidden_layers: int = 2,
                                            width_ratio: float = 5.0,
                                            convolutional_layers: int = 1,
                                            filter_size: int = 4,
                                            n_filters: int = 12,
                                            filter_stride: int = 3,
                                            n_encoded_sequence_features: int = 6,
                                            loss_fn=tf.losses.BinaryCrossentropy()):
        """
        Return a compiled neural network, similar to get_nn_model but also includes a convolutional network for
        encoding peptide sequences which feeds into the multilayer perceptron.

        :param learning_rate: Learning rate used by the optimizer (adam).
        :param dropout:  Dropout between each layer.
        :param hidden_layers: Number of hidden layers.
        :param width_ratio: Ratio of width of hidden layers to width of input layer.
        :param convolutional_layers: Number of convolutional layers.
        :param filter_size: Convolution filter size.
        :param n_filters: Number of filters.
        :param filter_stride: Filter stride.
        :param n_encoded_sequence_features: Number of nodes in the output of the convolutional network.
        :param loss_fn: The loss function to use.
        :return: A compiled keras.Model
        """
        optimizer = keras.optimizers.Adam(learning_rate=learning_rate)
        max_len = self.max_len if self.mhc_class == 'I' else self.max_len * 2
        model = get_model_with_peptide_encoding(ms_feature_length=self.feature_matrix.shape[1],
                                                dropout=dropout,
                                                hidden_layers_after_convolutions=hidden_layers,
                                                after_convolutions_width_ratio=width_ratio,
                                                convolutional_layers=convolutional_layers,
                                                filter_size=filter_size,
                                                n_filters=n_filters,
                                                filter_stride=filter_stride,
                                                n_encoded_sequence_features=n_encoded_sequence_features,
                                                max_pep_length=max_len
                                                )
        model.compile(optimizer=optimizer, loss=loss_fn)
        return model


    def run(self,
            model='MhcValidator',
            model_fit_function='fit',
            model_predict_function='predict',
            post_prediction_fn=lambda x: x,
            additional_training_data=None,
            return_prediction_data_and_model: bool = False,
            n_splits: int = 3,
            early_stopping_patience: int = 10,
            #q_value_subset: float = 1.0,
            #features_for_subset: Union[List[str], str] = 'all',
            #subset_threshold: int = 1,
            weight_by_inverse_peptide_counts: bool = False,
            visualize: bool = False,
            random_seed: int = None,
            clear_session: bool = True,
            alternate_labels=None,
            initial_model_weights: str = None,
            fit_model: bool = True,
            fig_pdf: Union[str, PathLike] = None,
            report_directory: Union[str, PathLike] = None,
            mhcflurry: bool = False,
            netmhcpan: bool = False,
            bigmhc: bool = False,
            mixmhc2pred: bool = False,
            peptdeep: bool = False,
            autort: bool = False,
            deeplc: bool = False,
            im2deep: bool = False,
            koina_predictors: List[str] = None,
            koina_server_url: str = 'koina.wilhelmlab.org:443',
            fine_tune: bool = False,
            sequence_encoding: bool = False,
            mzml_folder: PathLike = None,
            **kwargs):

        """
        Run the validation algorithm.

        :param model: The model to train on the target/decoy data. Can be a Python object with a fit and predict
        function, or string in {'NNValidator', 'MhcValidator'}. 'MhcValidator' causes predictions by MhcFlurry and
         NetMHCpan to be added to the feature matrix (if the respective arguments to `run` are set to True).
         NNValidator does not add the predictions.
        network supported by MhcValidator, while SEQUENCE_ENCODING will load the neural network which also performs
        convolutional peptide sequence encoding. Default is 'BASIC'.
        :param model_fit_function: The function which fits the model. Default is 'fit'.
        :param model_predict_function: The function used to make predictions with the fitted model. Default is 'predict'.
        :param post_prediction_fn: A function applied to the output of the predict function. Useful if the output is
        multidimensional, as all downstream processes expect a probability between 0 and 1.
        :param additional_training_data: Additional data for training the model. Only used if the model
        expects two inputs. If you are using the provided neural network which encodes peptide sequences, then you must
        pass self.X_encoded_peps.
        :param return_prediction_data_and_model: Whether to return predictions, q-values, etc in a dictionary. This
        data is available from the attributes of this MhcValidator instance after running, but it can be useful to
        return the data if you will be manipulating it downstream.
        :param n_splits: Number of splits used for training and validation/predicting (ala k-fold cross-validation).
        :param weight_by_inverse_peptide_counts: Whether to weight training by inverse peptide counts (i.e. number of
        times a sequence is identified in the data).
        :param visualize: Visualize the results.
        :param random_seed: Random seed used.
        :param clear_session: Clear the Tensorflow session before running.
        :param alternate_labels: Alternate labels to use for training. Possibly useful in an iterative variation of the
        algorithm.
        :param initial_model_weights: A file containing model weights to load before training using the models "load"
        function, if it has one.
        :param fit_model: Whether or not to fit the model. You would only set this to false if you were loading weights
        from an already-fitted model.
        :param fig_pdf: Filepath to save a PDF version of the training report.
        :param report_directory: Save all run information to a specified location. Includes: annotated input data,
        feature matrix, NetMHCpan and MHCFlurry predictions (if applicable), model weights, training report PDF.
        :param kwargs: Additional keyword arguments passed to model fit function.
        :return:
        """

        self.report_directory = str(report_directory)
        if not os.path.exists(report_directory):
            os.makedirs(report_directory, exist_ok=True)

        if clear_session:
            K.clear_session()

        if random_seed is None:
            random_seed = self.random_seed
        self._set_seed(random_seed)

        predictor_types = set()
        if autort or deeplc:
            predictor_types.add('RT')
        if im2deep:
            predictor_types.add('CCS')
        if peptdeep:
            predictor_types.add('RT') # TODO
            predictor_types.add('CCS')
            predictor_types.add('MS2')
        if koina_predictors:
            for predictor_name in koina_predictors:
                if predictor_name in KOINA_PREDICTORS.keys():
                    predictor_types.add(KOINA_PREDICTORS[predictor_name])
                    self.koina_predictors.append(predictor_name)
        if 'RT' in predictor_types and 'retentiontime' in self.raw_data:
            predictor_types.remove('RT')
            self.exp_rts = self.raw_data['retentiontime'].astype(float)
        if len(predictor_types) != 0:
            assert mzml_folder is not None, f'mzML folder must be provided for {predictor_types} scores'
            mzml_paths = Path(mzml_folder).rglob('*.mzML')
            mzml_map = {path.stem.replace('_uncalibrated', ''): str(path.expanduser().resolve()) for path in mzml_paths}
            mzml_name = self.filepath.name.replace('_edited.pin', '.pin').replace('.pin', '')
            assert mzml_name in mzml_map.keys(), f'mzML file not found: {mzml_folder}/{mzml_name}.mzML '
            mzml_path = mzml_map[mzml_name]
            if '_uncalibrated.' in mzml_path:
                self.exp_rts, self.exp_ims, self.exp_spectra = \
                    get_rt_ccs_ms2_from_msfragger_mzml(mzml_path, self.raw_data['ScanNr'].astype(int),
                                                       self.raw_data['ExpMass'].astype(float), self.charges)
            else:
                self.exp_rts, self.exp_ims, self.exp_spectra = \
                    get_rt_ccs_ms2_from_mzml(mzml_path, self.raw_data['ScanNr'].astype(int),
                                             self.raw_data['ExpMass'].astype(float), self.charges)

        self.fine_tune = fine_tune
        self.koina_server_url = koina_server_url
        if netmhcpan:
            self.add_netmhcpan_predictions()
        if mhcflurry and self.mhc_class == 'I':
            self.add_mhcflurry_predictions()
        if bigmhc and self.mhc_class == 'I':
            self.add_bigmhc_predictions()
        if mixmhc2pred and self.mhc_class == 'II':
            self.add_mixmhc2pred_predictions()
        if peptdeep:
            self.add_peptdeep_predictions()
        if autort:
            self.add_autort_predictions()
        if deeplc:
            self.add_deeplc_predictions()
        if im2deep:
            self.add_im2deep_predictions()
        if len(self.koina_predictors) > 0:
            self.add_koina_predictions()


        if model.lower() == 'mhcvalidator':
            if not sequence_encoding:
                model_args = {key: arg for key, arg in kwargs.items() if key in signature(self.get_nn_model).parameters}
                kwargs = {key: arg for key, arg in kwargs.items() if key not in model_args}
                model = self.get_nn_model(**model_args)
            else:
                model_args = {key: arg for key, arg in kwargs.items() if key in
                              signature(self.get_nn_model_with_sequence_encoding).parameters}
                kwargs = {key: arg for key, arg in kwargs.items() if key not in model_args}
                model = self.get_nn_model_with_sequence_encoding(**model_args)
                if self.encoded_peptides is None:
                    self.encode_peptide_sequences()
                additional_training_data = self.encoded_peptides

        if initial_model_weights is not None:
            model.load_weights(initial_model_weights)

        # check if the model is a Keras model, and if so check if number of epochs and batch size have been specified.
        # If they haven't set them to default values of 30 and 512, respectively. Otherwise things will go poorly.
        if isinstance(model, keras.Model):
            if 'epochs' not in kwargs.keys():
                print('`epochs` was not passed as a keyword argument. Setting it to default value of 30')
                kwargs['epochs'] = 30
            if 'batch_size' not in kwargs.keys():
                print('`batch_size` was not passed as a keyword argument. Setting it to default value of 512')
                kwargs['batch_size'] = 512

        # prepare data for training
        all_data = self.feature_matrix.copy(deep=True)

        if alternate_labels is None:
            labels = deepcopy(self.labels)
        else:
            labels = alternate_labels

        all_data = all_data.values
        peptides = self.peptides

        # we might make the splits better if our stratification takes feature q-values into account.
        # e.g. we calculate q-values for expect value and MHC predictions, and make sure we include good examples
        # from each allele.
        # stratification_labels = self.get_stratification_labels()
        '''skf = list(StratifiedKFold(n_splits=n_splits,
                                   random_state=random_seed,
                                   shuffle=True).split(all_data, labels))'''
        # skf = k_fold_split(peptides=peptides, k_folds=n_splits, random_state=random_seed)
        skf = k_fold_split_s(s=self.raw_data['ExpMass'].to_numpy(dtype=float), k_folds=n_splits, random_state=random_seed)
        # skf = k_fold_split_so(peptides=peptides, s=self.raw_data['ExpMass'].to_numpy(dtype=float), k_folds=n_splits, random_state=random_seed)

        predictions = np.zeros_like(labels, dtype=float)
        k_splits = np.zeros_like(labels, dtype=int)

        output = []
        history = []

        if isinstance(model, keras.Model):
            now = str(datetime.now()).replace(' ', '_').replace(':', '-')
            initial_model_weights = str(self.model_dir / f'mhcvalidator_initial_weights_{now}.keras')
            model.save(initial_model_weights)
        else:
            initial_model_weights = ''

        for k_fold, (train_index, predict_index) in enumerate(skf):
            print('-----------------------------------')
            print(f'Training on split {k_fold+1}')
            self._set_seed(random_seed)

            if isinstance(model, keras.Model):
                model.load_weights(initial_model_weights)
            feature_matrix = deepcopy(all_data)

            '''if q_value_subset < 1.:
                mask = self.get_qvalue_mask_from_features(X=feature_matrix[train_index],
                                                          y=labels[train_index],
                                                          cutoff=q_value_subset,
                                                          n=subset_threshold,
                                                          features_to_use=features_for_subset,
                                                          verbosity=1)
            else:
                mask = np.ones_like(labels[train_index], dtype=bool)'''
            mask = np.ones_like(labels[train_index], dtype=bool)  # just in case we implement the q-value subset again

            x_train = deepcopy(feature_matrix[train_index, :][mask])
            rnd_idx = RandomState(random_seed).choice(len(x_train), len(x_train), replace=False)
            x_train = x_train[rnd_idx]
            x_predict = deepcopy(feature_matrix[predict_index, :])
            input_scalar = NDStandardScaler()
            input_scalar = input_scalar.fit(x_train)
            x_train = input_scalar.transform(x_train)
            x_predict = input_scalar.transform(x_predict)
            feature_matrix = input_scalar.transform(feature_matrix)

            x = deepcopy(feature_matrix)
            x_train = deepcopy(x_train)
            x_predict = deepcopy(x_predict)
            train_labels = labels[train_index][mask][rnd_idx]
            predict_labels = labels[predict_index]
            print(f' Training split - {np.sum(train_labels == 1)} targets | {np.sum(train_labels == 0)} decoys')
            print(f' Prediction split - {np.sum(predict_labels == 1)} targets | {np.sum(predict_labels == 0)} decoys')

            if weight_by_inverse_peptide_counts:
                pep_counts = Counter(peptides[train_index][mask])
                weights = np.array([np.sqrt(1 / pep_counts[p]) for p in peptides[train_index][mask][rnd_idx]])
            else:
                weights = np.ones_like(labels[train_index][mask][rnd_idx])

            if additional_training_data is not None:
                additional_training_data = deepcopy(additional_training_data)
                x2_train = additional_training_data[train_index][mask][rnd_idx]
                x2_test = additional_training_data[predict_index]
                input_scalar2 = NDStandardScaler()
                input_scalar2 = input_scalar2.fit(x2_train)

                x2_train = input_scalar2.transform(x2_train)
                x2_test = input_scalar2.transform(x2_test)
                additional_training_data = input_scalar2.transform(additional_training_data)

                x_train = (x_train, x2_train)
                x_predict = (x_predict, x2_test)
                x = (x, additional_training_data)

            model_fit_parameters = eval(f'signature(model.{model_fit_function})').parameters
            if 'validation_data' in model_fit_parameters.keys():
                val_str = 'validation_data=(x_predict, predict_labels),'
            else:
                val_str = ''

            if 'sample_weight' in model_fit_parameters.keys():
                weight_str = 'sample_weight=weights,'
            else:
                weight_str = ''

            if isinstance(model, keras.Model):
                early_stopping = keras.callbacks.EarlyStopping(
                    monitor="val_loss",
                    patience=early_stopping_patience,
                    verbose=1,
                    mode="auto",
                    restore_best_weights=False)
                now = str(datetime.now()).replace(' ', '_').replace(':', '-')
                model_name = str(self.model_dir / f'mhcvalidator_k={k_fold+1}_{now}.keras')
                checkpoint = keras.callbacks.ModelCheckpoint(model_name,
                                                             monitor='val_loss', verbose=0,
                                                             save_best_only=True, mode='min')
                callbacks_str = 'callbacks=[early_stopping, checkpoint],'
            else:
                callbacks_str = ''
                model_name = ''

            # Train the model
            if fit_model:
                fit_history = eval(f"model.{model_fit_function}(x_train, train_labels, "
                                   f"{val_str} {weight_str} {callbacks_str} **kwargs)")
                if model_name != '':
                    model.load_weights(model_name)
                    if report_directory is not None:
                        if not os.path.exists(report_directory):
                            os.mkdir(report_directory)
                        model.save(Path(report_directory) / f'{Path(self.filename).stem}'
                                                            f'.mhcvalidator_model_k={k_fold+1}.keras')
            else:
                fit_history = None

            if fit_history is not None and hasattr(fit_history, 'history'):
                history.append(fit_history)

            predict_preds = post_prediction_fn(eval(
                f"model.{model_predict_function}(x_predict)")).flatten()  # all these predictions are assumed to be arrays. we flatten them because sometimes the have an extra dimension of size 1
            train_preds = post_prediction_fn(eval(f"model.{model_predict_function}(x_train)")).flatten()
            predict_qs = calculate_qs(predict_preds.flatten(), predict_labels)
            train_qs = calculate_qs(train_preds.flatten(), train_labels)
            preds = post_prediction_fn(eval(f"model.{model_predict_function}(x)")).flatten()
            qs = calculate_qs(preds.flatten(), labels)
            predictions[predict_index] = predict_preds
            k_splits[predict_index] = k_fold + 1
            assert np.all(predict_labels == self.labels[predict_index])

            train_roc = calculate_roc(train_qs, train_labels, qvalue_cutoff=0.05)
            val_roc = calculate_roc(predict_qs, predict_labels, qvalue_cutoff=0.05)
            roc = calculate_roc(qs, labels, qvalue_cutoff=0.05)

            pep_level_qs, _, pep_level_labels, peps, pep_counts = calculate_peptide_level_qs(predict_preds,
                                                                                             predict_labels,
                                                                                             self.peptides[predict_index])

            print(f' | PSMs in this split validated at 1% FDR: {np.sum((predict_qs <= 0.01) & (predict_labels == 1))}')
            print(f' | Peptides in this split validated at 1% FDR (peptide-level): '
                  f'{np.sum((pep_level_qs <= 0.01) & (pep_level_labels == 1))}')
            print('-----------------------------------')

            results = {'train_preds': train_preds, 'train_labels': train_labels, 'train_qs': train_qs,
                           'train_roc': train_roc, 'predict_preds': predict_preds, 'predict_labels': predict_labels,
                           'predict_qs': predict_qs,
                           'predict_roc': val_roc, 'preds': preds, 'labels': labels, 'qs': qs, 'roc': roc, 'model': model,
                       'train_index': train_index, 'predict_index': predict_index}
            output.append(results)

        self.predictions = np.empty(len(labels), dtype=float)
        self.qs = np.empty(len(labels), dtype=float)

        self.predictions = predictions
        self.qs = calculate_qs(predictions, labels)
        self.roc = calculate_roc(self.qs, self.labels)

        pep_level_qs, _, pep_level_labels, pep_level_peps, pep_counts = calculate_peptide_level_qs(self.predictions,
                                                                                                   self.labels,
                                                                                                   self.peptides)

        print('===================================')
        print('Validation results')
        print(f' | PSMs validated at 1% FDR: {np.sum((self.qs <= 0.01) & (self.labels == 1))}')
        print(f' | Peptides validated at 1% FDR (peptide-level): '
              f'{np.sum((pep_level_qs <= 0.01) & (pep_level_labels == 1))}')
        print('===================================')

        fig = plt.figure(constrained_layout=True, figsize=(10, 10))
        fig.suptitle(self.filename, fontsize=16)
        gs = GridSpec(2, 2, figure=fig)

        # train = fig.add_subplot(gs[:4, 0])
        # val = fig.add_subplot(gs[4:8, 0])
        final: plt.Axes = fig.add_subplot(gs[0, 0])

        if len(history) > 0:  # the model returned a fit history we can use here
            dist: plt.Axes = fig.add_subplot(gs[0, 1])
            loss: plt.Axes = fig.add_subplot(gs[1, 1])
            #train_split = fig.add_subplot(gs[8:10, 1])
            #val_split = fig.add_subplot(gs[10:, 1])
        else:
            loss: plt.Axes = None
            dist: plt.Axes = fig.add_subplot(gs[0, 1])
            #train_split = fig.add_subplot(gs[6:9, 1])
            #val_split = fig.add_subplot(gs[9:, 1])

        colormap = get_cmap("tab10")

        # self._visualize_splits(skf, split='train', ax=train_split)
        # self._visualize_splits(skf, split='val', ax=val_split)
        # train_split.set_title('K-fold splits')
        # train_split.set_ylabel('Training')
        # val_split.set_ylabel('Validation')
        # val_split.set_xlabel('Scan number')

        if loss:
            min_x = []
            min_y = []
            for i, h in enumerate(history):
                loss.plot(range(1, len(h.history['val_loss']) + 1),
                          h.history['val_loss'], c=colormap(i), marker=None, label=f'split {i+1}')
                min_y.append(np.min(h.history['val_loss']))
                min_x.append(np.argmin(h.history['val_loss']) + 1)
            loss.plot(min_x, min_y, ls='none', marker='x', ms='12', c='k', label='best models')
            loss.set_title('Validation loss')
            loss.set_xlabel('Epoch')
            loss.set_ylabel('Loss')
            loss.legend()

        # for i, r in enumerate(output):
        #     train.plot(*r['train_roc'], c=colormap(i), ms='3', ls='none', marker='.', label=f'split {i+1}', alpha=0.6)
        #     val.plot(*r['predict_roc'], c=colormap(i), ms='3', ls='none', marker='.', label=f'split {i+1}', alpha=0.6)
        final.plot(*self.roc, c=colormap(0), ms='3', ls='none', marker='.', alpha=0.6)
        n_psms_at_1percent = np.sum((self.qs <= 0.01) & (self.labels == 1))
        final.vlines(0.01, 0, n_psms_at_1percent, ls='--', lw=1, color='k', alpha=0.7)
        final.hlines(n_psms_at_1percent, 0, 0.01, ls='--', lw=1, color='k', alpha=0.7)

        _, bins, _ = dist.hist(self.predictions[self.labels == 1], label='Target', bins=30, alpha=0.5, color='g')
        dist.hist(self.predictions[self.labels == 0], label='Decoy', bins=bins, alpha=0.5, zorder=100, color='r')

        # train.set_xlim((0, 0.05))
        # val.set_xlim((0, 0.05))
        final.set_xlim((0, 0.05))

        # train.set_title('Training data')
        # train.set_xlabel('q-value')
        # train.set_ylabel('PSMs')
        # train.set_ylim((0, train.get_ylim()[1]))

        # val.set_title('Validation data')
        # val.set_xlabel('q-value')
        # val.set_ylabel('PSMs')
        # val.set_ylim((0, val.get_ylim()[1]))

        final.set_title('Final q-values')
        final.set_xlabel('q-value')
        final.set_ylabel('PSMs')
        final.set_ylim((0, final.get_ylim()[1]))

        dist.set_title('Prediction distributions')
        dist.set_xlabel('Target probability')
        dist.set_ylabel('PSMs')

        # train.legend(markerscale=3)
        # val.legend(markerscale=3)
        dist.legend()

        plt.tight_layout()
        if fig_pdf is not None:
            pdf = plt_pdf.PdfPages(str(fig_pdf), keep_empty=False)
            pdf.savefig(fig)
            pdf.close()
        if report_directory is not None:
            pdf_file = Path(report_directory) / f'{Path(self.filename).stem}.MhcValidator_training_report.pdf'
            pdf = plt_pdf.PdfPages(str(pdf_file), keep_empty=False)
            pdf.savefig(fig)
            pdf.close()
        if visualize:
            fig.show()
        plt.close(fig)

        # make peptide-level q-value lookup
        pep_q_lookup = {pep: q for pep, q in zip(pep_level_peps, pep_level_qs)}

        self.raw_data['mhcv_peptide'] = self.peptides
        self.raw_data['mhcv_prob'] = self.predictions
        self.raw_data['mhcv_label'] = self.labels
        self.raw_data['mhcv_q-value'] = self.qs
        self.raw_data['mhcv_pep-level_q-value'] = np.array([pep_q_lookup[p] for p in self.peptides])
        self.raw_data['mhcv_k-fold_split'] = k_splits

        self.annotated_data = self.raw_data.copy(deep=True)

        if report_directory is not None:
            self.annotated_data.to_csv(Path(report_directory) /
                                       f'{Path(self.filename).stem}.MhcValidator_annotated.tsv',
                                       index=False, sep='\t')
            features_all = self.raw_data.join(self.feature_matrix, how='left', rsuffix='_right')
            features_all = features_all[[col for col in features_all.columns if '_right' not in col]]
            features_all.to_csv(Path(report_directory) / f'{Path(self.filename).stem}.features.tsv',
                                       index=False, sep='\t')
            if self.mhcflurry_predictions is not None:
                self.mhcflurry_predictions.to_csv(Path(report_directory) /
                                                  f'{Path(self.filename).stem}.MhcFlurry_Predictions.tsv',
                                                  index=False, sep='\t')
            if self.netmhcpan_predictions is not None:
                self.netmhcpan_predictions.to_csv(Path(report_directory) /
                                                  f'{Path(self.filename).stem}.NetMHCpan_Predictions.tsv',
                                                  index=False, sep='\t')

        if return_prediction_data_and_model:
            return output, {'predictions': deepcopy(self.predictions),
                            'qs': deepcopy(self.qs),
                            'roc': deepcopy(self.roc)}



    def get_peptide_list_at_fdr(self, fdr: float, label: int = 1, peptide_level: bool = False):
        if peptide_level:
            qs, _, labels, peps, _ = calculate_peptide_level_qs(self.predictions, self.labels,
                                                                self.peptides)
            return peps[(qs <= fdr) & (labels == label)]
        else:
            return self.peptides[(self.qs <= fdr) & (self.labels == label)]
