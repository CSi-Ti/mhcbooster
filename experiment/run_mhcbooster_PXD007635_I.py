
import sys
from pathlib import Path

from mhcnames import normalize_allele_name

sys.path.append('/home/rw762/workspace/mhc-booster')
from mhcvalidator.validator import MhcValidator

allele_map_I = {
    'OvCa9' :   ['A*02:01', 'A*03:01', 'B*07:02', 'B*40:02', 'C*07:02', 'C*12:01'],
    'OvCa10' :  ['A*02:01', 'A*11:01', 'B*44:05', 'B*51:01', 'C*02:02', 'C*15:02'],
    'OvCa12' :  ['A*24:02', 'A*31:01', 'B*35:03', 'B*49:01', 'C*07:01', 'C*12:03'],
    'OvCa13' :  ['A*02', 'B*35', 'B*40', 'C*03', 'C*04'],   ###
    'OvCa15' :  ['A*11:01', 'A*24:02', 'B*07:02', 'B*55:01', 'C*03:03', 'C*07:02'],
    'OvCa16' :  ['A*02', 'B*40', 'B*4402', 'C*03', 'C*05'],   ###
    'OvCa23' :  ['A*01', 'A*03', 'B*08', 'B*35', 'C*04', 'C*07'],   ###
    'OvCa28' :  ['A*01:01', 'A*02:01', 'B*27:05', 'B*52:01', 'C*01:02', 'C*02:02'],
    'OvCa39' :  ['A*25:01', 'A*31:01', 'B*07:02', 'B*18:01', 'C*12:03', 'C*07:02'],
    'OvCa41' :  ['A*02', 'A*2402', 'B*18', 'B*51', 'C*02', 'C*12'],   ###
    'OvCa43' :  ['A*02', 'A*32', 'B*18', 'B*35', 'C*04', 'C*07'],   ###
    'OvCa45' :  ['A*01', 'A*23', 'B*08', 'B*44', 'C*04', 'C*07'],   ###
    'OvCa48' :  ['A*02:01', 'A*25:01', 'B*15:01', 'B*41:02', 'C*03:04', 'C*17:01'],
    'OvCa53' :  ['A*02', 'A*03', 'B*27', 'B*35', 'C*02', 'C*04'],   ###
    'OvCa54' :  ['A*0201', 'A*11:01', 'B*35:01', 'B*35:03', 'C*04:01', 'C*12:03'],
    'OvCa57' :  ['A*25', 'A*32', 'B*15', 'B*18', 'C*03', 'C*12'],   ###
    'OvCa58' :  ['A*02', 'A*03', 'B*35', 'C*03', 'C*04'],   ###
    'OvCa59' :  ['A*03', 'A*30', 'B*13', 'C*06'],   ###
    'OvCa60' :  ['A*24:02', 'A*25:01', 'B*13:02', 'B*18:01', 'C*12:03', 'C*06:02'],
    'OvCa64' :  ['A*01', 'A*25', 'B*08', 'C*07'],   ###
    'OvCa65' :  ['A*01', 'A*2402', 'B*15', 'B*35', 'C*04', 'C*14'],   ###
    'OvCa66' :  ['A*11:01', 'A*29:02', 'B*18:01', 'B*44:03', 'C*05:01', 'C*16:01'],
    'OvCa68' :  ['A*02:01', 'A*01:01', 'B*44:02', 'B*37:01', 'C*06:02', 'C*05:01'],
    'OvCa69' :  [],
    'OvCa70' :  ['A*01', 'A*02', 'B*07', 'C*07'],   ###
    'OvCa72' :  ['A*03:01', 'A*01:01', 'B*08:01', 'B*07:02', 'C*07:02', 'C*07:01'],
    'OvCa73' :  ['A*01:01', 'B*08:01', 'C*07:01' ],
    'OvCa74' :  ['A*02:01', 'B*18:01', 'B*51:01', 'C*07:02', 'C*15:02'],
    'OvCa79' :  ['A*01:01', 'A*31:01', 'B*08:01', 'B*51:01', 'C*07:01', 'C*15:02'],
    'OvCa80' :  ['A*25:01', 'A*32:01', 'B*18:01', 'B*39:01', 'C*12:03'],
    'OvCa81' :  ['A*02:01', 'B*45:01', 'B*56:01', 'C*07:02', 'C*01:02'],
    'OvCa82' :  ['A*01:01', 'A*03:01', 'B*08:01', 'B*38:01', 'C*07:01', 'C*12:03'],
    'OvCa83' :  ['A*02', 'A*11', 'B*51', 'B*55', 'C*03', 'C*15'],   ###
    'OvCa84' :  ['A*02:01', 'B*07:02', 'B*44:02', 'C*07:02', 'C*05:01'],
    'OvCa99' :  ['A*02:01', 'A*24:02', 'B*13:02', 'B*40:01', 'C*03:04', 'C*06:02'],
    'OvCa100' : ['A*02:01', 'B*07:02', 'B*41:02', 'C*07:02', 'C*17:01'],
    'OvCa103' : ['A*02:01', 'A*24:02', 'B*27:02', 'B*27:05', 'C*02:02'],
    'OvCa104' : ['A*03:01', 'B*07:02', 'B*35:08', 'C*04:01', 'C*07:02'],
    'OvCa105' : ['A*26:01', 'A*68:01', 'B*18:01', 'B*55:01', 'C*03:03',' C*07:01'],
    'OvCa107' : ['A*02:05', 'A*25:01', 'B*07:02', 'B*14:02', 'C*07:02',' C*08:02'],
    'OvCa109' : ['A*02:01', 'A*23:01', 'B*40:01', 'B*49:01', 'C*07:01',' C*03:04'],
    'OvCa111' : ['A*01:01', 'A*25:01', 'B*08:01', 'B*44:02', 'C*05:01',' C*07:01'],
    'OvCa114' : ['A*29:02', 'B*44:03', 'C*16:01']
}

allele_map_II = { }

allele_map = allele_map_I

pin_file_dir = '/mnt/d/workspace/mhc-validator-2/data/PXD007635/HLA-I'
mzml_folder = Path('/mnt/d/workspace/mhc-validator-2/data/PXD007635/HLA-I')
output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/PXD007635/HLA-I/test')
output_folder.mkdir(parents=True, exist_ok=True)

avail_alleles = [line.split(' ')[0] for line in open('/mnt/d/workspace/mhc-validator-2/third_party/netMHCpan-4.1/Linux_x86_64/data/MHC_pseudo.dat')]

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
    if len(alleles) == 0 or 'CLL' in file_name:
        use_mhc_scores = False
        print('Warning: no allele specified for {}'.format(file_name))

    print(file_name)

    validator = MhcValidator(max_threads=20)  # Open a MHCvalidator instance, a new one has to be opened for each .pin file
    validator.set_mhc_params(alleles=alleles, mhc_class='I', min_pep_len=8)  # Load the alleles you specified above
    validator.load_data(pin, filetype='pin')  # Load the pin file
    validator.run(sequence_encoding=True, mhcflurry=use_mhc_scores, bigmhc=use_mhc_scores,
                  netmhcpan=use_mhc_scores, mixmhc2pred=use_mhc_scores,
                  autort=True, deeplc=False,
                  im2deep=False,
                  peptdeep=False,
                  # koina_predictors=rt_all_models + ms2_all_models + ccs_top_models,
                  koina_predictors=['Deeplc_hela_hf', 'ms2pip_timsTOF2024', 'ms2pip_timsTOF2023'],
                  # koina_predictors=['Prosit_2024_irt_cit', 'Prosit_2019_irt', 'ms2pip_timsTOF2024', 'Prosit_2023_intensity_timsTOF', 'AlphaPeptDeep_ccs_generic'],
                  # koina_predictors=['Prosit_2024_irt_cit', 'Prosit_2019_irt'],
                  fine_tune=False,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  report_directory=output_folder / f'{file_name}')
