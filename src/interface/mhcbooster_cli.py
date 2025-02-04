import os
import re
import sys
import argparse

from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src import __version__
from src.main_mhcbooster import MhcValidator


description = f"""
MhcValidator v{__version__} (https://github.com/CaronLab/???)
Copyright 2024 Ruimin Wang under GNU General Public License v3.0

MhcValidator is a tool for validating peptide-spectrum matches from mass spectrometry database searches. It is 
intended for use with data from immunopeptidomics experiments, though it can be use for most types of 
proteomics experiments as well.
"""

parser = argparse.ArgumentParser(description=description)

general = parser.add_argument_group('general parameters')
general.add_argument('-i',
                     '--input',
                     required=True,
                     nargs='+',
                     type=str,
                     help='Input file(s) for MhcValidator. Must be comma- or tab-separated files or pepXML. Note that '
                          'MhcValidator has only been thoroughly tested using PIN files as input '
                          '(Percolator input files). You can pass multiple files as a space-separated list. If you '
                          'pass a generic tabular file, it must contain a column titled "Peptide" or "peptide" which '
                          'contains the peptide sequences. For generic tabular files, you should also use the '
                          '--prot_column, --decoy_tag, --tag_is_prefix arguments so '
                          'MhcValidator can figure out which PSMs are targets and which are decoys.')

general.add_argument('-m',
                     '--mzml_dir',
                     type=str,
                     help='Directory for mzML files, which are mandatory for MS2 & CCS score calculation. '
                          'It is also mandatory if the PSM files do not contain retention time column. '
                          'The _uncalibrated.mzML files from MSFragger are recommended.')

general.add_argument('-o',
                     '--output_dir',
                     type=str,
                     help='Output directory for MhcValidator. If not indicated, the input directory will be used.')

general.add_argument('--pep_column',
                     type=str,
                     help='The header of the column containing peptide sequences. Generally not required unless '
                          'the input is a generic text file (e.g. a CSV export from a search engine).')

general.add_argument('--prot_column',
                     type=str,
                     help='The header of the column containing protein identifications. Used '
                          'for inferring which PSMs are targets and which are decoys. Generally not required unless '
                          'the input is a generic text file (e.g. a CSV export from a search engine).')

general.add_argument('--decoy_tag',
                     type=str,
                     help='The tag indicating decoy hits in the protein column, e.g. rev_ or decoy_ are common. Used '
                          'for inferring which PSMs are targets and which are decoys. Usually not required for '
                          'PIN files.')

general.add_argument('--tag_is_prefix',
                     type=bool,
                     default=True,
                     help='Whether the decoy tag is a prefix or not. If not, it is assumed to be a suffix. Used '
                          'for inferring which PSMs are targets and which are decoys. Usually not required for '
                          'PIN files.')

general.add_argument('--delimiter',
                     type=str,
                     default='\t',
                     help='The delimiter of the file, if it is tabular data.')

general.add_argument('--min_pep_len',
                     type=int,
                     default=8,
                     help='The minimum peptide length to consider.')

general.add_argument('--max_pep_len',
                     type=int,
                     default=30,
                     help='The maximum peptide length to consider.')

general.add_argument('-n',
                     '--n_processes',
                     type=int,
                     default=0,
                     help='The number of threads to be used concurrently when running NetMHCpan. Uses all available '
                          'CPUs if < 1.')

mhc_params = parser.add_argument_group('Rescoring parameters', 'MHC/RT/MS2/CCS prediction parameters.')

mhc_params.add_argument('-a',
                        '--alleles',
                        nargs='+',
                        type=str,
                        help='MHC allele(s) of the sample of interest. If there is more than one, pass them as a space-'
                             'separated list. Not required if you are not running APP predictors.')
mhc_params.add_argument('-c',
                        '--mhc_class',
                        type=str,
                        choices=('I', 'II'),
                        help='The class of MHC allele(s) and predictors.')

mhc_params.add_argument('-app',
                        '--app_predictors',
                        nargs='+',
                        type=str,
                        choices=('NetMHCpan', 'MHCflurry', 'BigMHC', 'NetMHCIIpan', 'MixMHC2pred'),
                        help='The APP score predictors you want to be considered by the discriminant function.')

mhc_params.add_argument('-rt',
                        '--rt_predictors',
                        nargs='+',
                        type=str,
                        choices=('AutoRT', 'DeepLC', 'Deeplc_hela_hf', 'AlphaPeptDeep_rt_generic', 'Prosit_2019_irt',
                                 'Prosit_2024_irt_cit', 'Prosit_2020_irt_TMT', 'Chronologer_RT'),
                        help='The RT score predictors you want to be considered by the discriminant function.')

mhc_params.add_argument('-ms2',
                        '--ms2_predictors',
                        nargs='+',
                        type=str,
                        choices=('AlphaPeptDeep_ms2_generic', 'Prosit_2019_intensity', 'Prosit_2024_intensity_cit',
                                 'Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_CID',
                                 'Prosit_2020_intensity_HCD', 'Prosit_2020_intensity_TMT',
                                 'ms2pip_HCD2021', 'ms2pip_timsTOF2023', 'ms2pip_iTRAQphospho', 'ms2pip_Immuno_HCD',
                                 'ms2pip_TTOF5600', 'ms2pip_timsTOF2024', 'ms2pip_CID_TMT'),
                        help='The MS2 score predictors you want to be considered by the discriminant function.')

mhc_params.add_argument('-ccs',
                        '--ccs_predictors',
                        nargs='+',
                        type=str,
                        choices=('IM2Deep', 'AlphaPeptDeep_ccs_generic'),
                        help='The RT score predictors you want to be considered by the discriminant function.')

mhc_params.add_argument('--koina_server_url',
                        type=str,
                        default='koina.wilhelmlab.org:443',
                        help='The URL of Koina server for RT, MS2 and CCS prediction. Default server is koina.wilhelmlab.org:443')


mhc_params.add_argument('--fine_tune',
                        action='store_true',
                        help='Fine-tune the models before prediction. Supported models: [AutoRT, DeepLC, IM2Deep]')

mhc_params.add_argument('--auto_pred',
                        action='store_true',
                        help='If you cannot decide which predictors to use, try use --auto-pred. '
                             'This function will predict the best combination of predictors. '
                             'Predictors will be fully override by the predicted combination.')

training = parser.add_argument_group('training parameters', 'Related to the training of the artificial neural network.')

training.add_argument('-v',
                      '--verbose_training',
                      type=int,
                      default=0,
                      help='The verbosity level of tensorflow during training. Should be one of {0, 1, 2}.')

training.add_argument('-k',
                      '--k_folds',
                      type=int,
                      default=5,
                      help='The number of splits used in training and predictions, as in K-fold cross-validation.')

training.add_argument('-s',
                      '--encode_peptide_sequences',
                      action='store_true',
                      help='Encode peptide sequences as features for the training algorithm.')

def run():
    args = parser.parse_args()

    input_files = args.input
    if len(input_files) == 1 and os.path.isdir(input_files[0]):
        input_files = Path(input_files[0]).rglob('*.pin')

    alleles = []
    allele_map = {}
    if args.alleles is not None:
        if len(args.alleles) == 1 and os.path.exists(args.alleles[0]):
            for line in open(args.alleles[0]):
                line_split = re.split(r'[\t,]', line)
                allele_map[line_split[0].strip()] = [allele.strip() for allele in line_split[1].split(';')]
        else:
            alleles = args.alleles

    for input_file in input_files:
        print(f'Processing: {input_file}')

        if args.output_dir is None:
            args.output_dir = Path(input_file).parent

        run_alleles = alleles.copy()
        if len(run_alleles) == 0 and len(allele_map) != 0:
            for keyword in allele_map.keys():
                if keyword in input_file.stem:
                    run_alleles = allele_map[keyword]
                    break

        v = MhcValidator(max_threads=args.n_processes)
        v.set_mhc_params(alleles=run_alleles, mhc_class=args.mhc_class, min_pep_len=args.min_pep_len, max_pep_len=args.max_pep_len)
        v.load_data(input_file,
                    peptide_column=args.pep_column,
                    protein_column=args.prot_column,
                    decoy_tag=args.decoy_tag,
                    tag_is_prefix=args.tag_is_prefix,
                    file_delimiter=args.delimiter)

        v.run(sequence_encoding=args.encode_peptide_sequences,
              app_predictors=args.app_predictors,
              rt_predictors=args.rt_predictors,
              ms2_predictors=args.ms2_predictors,
              ccs_predictors=args.ccs_predictors,
              auto_predict_predictor=args.auto_pred,
              fine_tune=args.fine_tune,
              mzml_folder=args.mzml_dir,
              report_directory=Path(args.output_dir) / f'{Path(input_file).stem}_MhcValidator',
              n_splits=args.k_folds,
              visualize=False,
              verbose=args.verbose_training)

        if args.auto_pred:
            args.rt_predictors = v.rt_predictors
            args.ms2_predictors = v.ms2_predictors
            args.ccs_predictors = v.ccs_predictors
            args.app_predictors = v.app_predictors
            args.auto_pred = False



if __name__ == '__main__':
    run()
