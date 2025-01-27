
import sys
from pathlib import Path

sys.path.append('/home/rw762/workspace/mhc-booster')
from mhcvalidator.validator import MhcValidator

allele_map_I = {
    'HLA': ['A*0101', 'A*0202', 'B*5701', 'B*4403', 'C*0602', 'C*1602'],
}

allele_map_II = {
    'HLA': ['DQA1*03:02', 'DQA1*02:01', 'DQB1*03:03', 'DRB1*04:05', 'DRB1*07:01']
}

allele_map = allele_map_I

pin_file_dir = '/mnt/d/workspace/mhc-validator-2/data/MSV000091456/A375_lowInput_IP/HLA-I'
mzml_folder = Path('/mnt/d/workspace/mhc-validator-2/data/MSV000091456/A375_lowInput_IP/HLA-I')
output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/MSV000091456/A375_lowInput_IP/HLA-I/test')
output_folder.mkdir(parents=True, exist_ok=True)

pin_files = list(Path(pin_file_dir).rglob('*.pin'))

rt_all_models = ['Deeplc_hela_hf', 'AlphaPeptDeep_rt_generic', 'Prosit_2019_irt', 'Prosit_2024_irt_cit']
ms2_all_models = ['Prosit_2019_intensity', 'Prosit_2024_intensity_cit', 'Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_CID', 'Prosit_2020_intensity_HCD',
                  'ms2pip_HCD2021', 'ms2pip_timsTOF2023', 'ms2pip_Immuno_HCD', 'ms2pip_timsTOF2024']
ccs_top_models = ['AlphaPeptDeep_ccs_generic']

for i in range(len(pin_files)):
    print(f'TracerID: {i}')
    pin = pin_files[i]
    file_name = pin.stem
    if 'edited' in file_name:
        continue
    alleles = []
    for donor in allele_map.keys():
        if donor in file_name:
            alleles = allele_map[donor]
            break
    use_mhc_scores = True

    print(file_name)

    validator = MhcValidator(max_threads=20)  # Open a MHCvalidator instance, a new one has to be opened for each .pin file
    validator.set_mhc_params(alleles=alleles, mhc_class='I', min_pep_len=8)  # Load the alleles you specified above
    validator.load_data(pin, filetype='pin')  # Load the pin file
    validator.run(sequence_encoding=True, mhcflurry=use_mhc_scores, bigmhc=use_mhc_scores,
                  netmhcpan=use_mhc_scores, mixmhc2pred=use_mhc_scores,
                  autort=False, deeplc=False,
                  im2deep=True,
                  peptdeep=False,
                  # koina_predictors=rt_all_models + ms2_all_models + ccs_top_models,
                  koina_predictors=['Deeplc_hela_hf', 'Prosit_2019_irt', 'AlphaPeptDeep_rt_generic',
                                    'ms2pip_timsTOF2024', 'ms2pip_timsTOF2023', 'AlphaPeptDeep_ccs_generic'],
                  fine_tune=False,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  report_directory=output_folder / f'{file_name}')
