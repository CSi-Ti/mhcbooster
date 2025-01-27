import os
import sys
import argparse

from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mhcvalidator import __version__
from mhcvalidator.validator import MhcValidator


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
                             'separated list. Not required if you are not running MhcFlurry or NetMHCpan.')

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
                                 'Prosit_2024_irt_cit', 'Chronologer_RT'),
                        help='The RT score predictors you want to be considered by the discriminant function.')

mhc_params.add_argument('-ms2',
                        '--ms2_predictors',
                        nargs='+',
                        type=str,
                        choices=('Prosit_2019_intensity', 'Prosit_2024_intensity_cit', 'Prosit_2023_intensity_timsTOF',
                                 'Prosit_2020_intensity_CID', 'Prosit_2020_intensity_HCD',
                                 'ms2pip_HCD2021', 'ms2pip_timsTOF2023', 'ms2pip_Immuno_HCD', 'ms2pip_timsTOF2024'),
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

    for input_file in input_files:
        if args.output_dir is None:
            args.output_dir = Path(input_file).parent


        use_netmhcpan, use_mhcflurry, use_bigmhc, use_netmhcIIpan, use_mixmhc2pred = False, False, False, False, False
        if args.app_predictors is not None:
            use_netmhcpan = 'NetMHCpan' in args.app_predictors
            use_mhcflurry = 'MHCflurry' in args.app_predictors
            use_bigmhc = 'BigMHC' in args.app_predictors
            use_netmhcIIpan = 'NetMHCIIpan' in args.app_predictors
            use_mixmhc2pred = 'MixMHC2pred' in args.app_predictors
        mhc_class = 'I'
        if use_netmhcIIpan or use_mixmhc2pred:
            mhc_class = 'II'

        use_autort, use_deeplc = False, False
        koina_rt_predictors = []
        if args.rt_predictors is not None:
            use_autort = 'AutoRT' in args.rt_predictors
            use_deeplc = 'DeepLC' in args.rt_predictors
            koina_rt_predictors = [m for m in args.rt_predictors if m != 'AutoRT' and m != 'DeepLC']

        use_im2deep = False
        koina_ccs_predictors = []
        if args.ccs_predictors is not None:
            use_im2deep = 'IM2Deep' in args.ccs_predictors and args.fine_tune
            koina_ccs_predictors = args.ccs_predictors
        if use_im2deep:
            koina_ccs_predictors = [m for m in koina_ccs_predictors if m != 'IM2Deep']

        koina_ms2_predictors = []
        if args.ms2_predictors is not None:
            koina_ms2_predictors = args.ms2_predictors

        v = MhcValidator(max_threads=args.n_processes)
        v.set_mhc_params(alleles=args.alleles, mhc_class=mhc_class, min_pep_len=args.min_pep_len, max_pep_len=args.max_pep_len)
        v.load_data(input_file,
                    peptide_column=args.pep_column,
                    protein_column=args.prot_column,
                    decoy_tag=args.decoy_tag,
                    tag_is_prefix=args.tag_is_prefix,
                    file_delimiter=args.delimiter)

        v.run(sequence_encoding=args.encode_peptide_sequences,
              netmhcpan=use_netmhcpan or use_netmhcIIpan, mhcflurry=use_mhcflurry, bigmhc=use_bigmhc, mixmhc2pred=use_mixmhc2pred,
              autort=use_autort, deeplc=use_deeplc,
              im2deep=use_im2deep,
              koina_predictors=koina_rt_predictors + koina_ms2_predictors + koina_ccs_predictors,
              fine_tune=args.fine_tune,
              mzml_folder=args.mzml_dir,
              report_directory=Path(args.output_dir) / f'{Path(input_file).stem}_MhcValidator',
              n_splits=args.k_folds,
              visualize=False,
              verbose=args.verbose_training)


if __name__ == '__main__':
    run()
