
import sys
from pathlib import Path

from mhcnames import normalize_allele_name

sys.path.append('/home/rw762/workspace/mhc-booster')
from mhcvalidator.validator import MhcValidator

allele_map_I = {
    'UPN001' :  ['A*02:01',	'A*11:01',	'B*39:01',	'B*40:01',	'C*03:04',	'C*12:03'],
    'UPN002' :  ['A*02:01',	'B*35:01',	'B*39:01',	'C*04:01',	'C*12:03'],
    'UPN003' :  ['A*25:01',	'A*26:01',	'B*18:01',	'B*38:01',	'C*12:03'],
    'UPN004' :  ['A*01:01',	'A*24:02',	'B*08:01',	'B*27:05',	'C*02:02',	'C*07:02'],
    'UPN005' :  ['A*02:01',	'A*03:01',	'B*18:01',	'B*35:01',	'C*04:01',	'C*05:01'],
    'UPN006' :  ['A*03:01',	'A*30:01',	'B*07:02',	'B*13:02',	'C*06:02',	'C*07:02'],
    'UPN007' :  ['A*02:01',	'A*03:01',	'B*07:02',	'B*55:01',	'C*01:02',	'C*07:02'],
    'UPN008' :  ['A*01:01',	'A*02:01',	'B*27:05',	'B*37:01',	'C*02:02',	'C*06:02'],
    'UPN009' :  ['A*01:01',	'A*02:01',	'B*27:02',	'B*37:01',	'C*02:02',	'C*06:02'],
    'UPN010' :  ['A*23:01',	'B*49:01',	'C*03:04',	'C*07:02'],
    'UPN011' :  ['A*01:01',	'A*68:01',	'B*08:01',	'B*44:02',	'C*07:01',	'C*07:04'],
    'UPN012' :  ['A*02:01',	'A*03:01',	'B*40:01',	'C*03:04'],
    'UPN013' :  ['A*01:01',	'A*03:01',	'B*08:01',	'B*51:01',	'C*01:02',	'C*07:01'],
    'UPN014' :  ['A*02:01',	'A*03:01',	'B*07:02',	'B*44:02',	'C*05:01',	'C*07:02'],
    'UPN015' :  ['A*02:05',	'A*24:02',	'B*35:01',	'B*50:01',	'C*04:01',	'C*06:02'],
    'UPN016' :  ['A*11:01',	'B*15:01',	'B*52:01',	'C*04:01',	'C*12:02'],
    'UPN017' :  ['A*01:01',	'A*02:01',	'B*08:01',	'B*13:02',	'C*06:02',	'C*07:02'],
    'UPN018' :  ['A*32:01',	'A*68:01',	'B*07:02',	'B*27:05',	'C*07:01',	'C*07:02'],
    'UPN019' :  ['A*01:01',	'A*26:01',	'B*07:02',	'B*37:01',	'C*02:02',	'C*07:02'],
    'UPN020' :  ['A*01:01',	'A*32:01',	'B*07:02',	'B*44:02',	'C*07:01',	'C*07:02'],
    'UPN021' :  ['A*02:01',	'A*24:02',	'B*51:01',	'B*57:01',	'C*14:02'],
    'UPN022' :  ['A*02:01',	'A*03:01',	'B*40:01',	'B*44:02',	'C*02:02',	'C*03:04'],
    'UPN023' :  ['A*02:01',	'A*03:01',	'B*35:01',	'C*04:01'],
    'UPN024' :  ['A*01:01',	'A*02:01',	'B*08:01',	'B*51:01',	'C*02:02',  'C*07:02'],
    'UPN025' :  ['A*02:01',	'A*24:02',	'B*51:01',	'B*57:01',	'C*01:02',	'C*06:02'],
    'UPN026' :  ['A*01:01',	'A*03:01',	'B*07:02',	'B*44:02',	'C*05:01',	'C*07:02'],
    'UPN027' :  ['A*01:01',	'A*02:01',	'B*08:01',	'B*27:02',	'C*02:02',	'C*07:01'],
    'UPN028' :  ['A*02:01',	'B*15:01',	'B*56:01',	'C*01:02',	'C*03:04'],
    'UPN029' :  ['A*24:02',	'A*26:01',	'B*27:05',	'B*39:01',	'C*01:02',	'C*07:02'],
    'UPN030' :  ['A*02:01',	'B*07:02',	'B*18:01',	'C*03:04',	'C*06:02'],
    'UPN031' :  ['A*02:01',	'A*11:01',	'B*35:01',	'B*40:01',	'C*03:04',	'C*04:01'],
    'UPN032' :  ['A*01:01',	'A*24:02',	'B*07:02',	'B*08:01',	'C*07:02'],
    'UPN033' :  ['A*02:01',	'A*68:01',	'B*38:01',	'B*51:01',	'C*12:03',	'C*14:02'],
    'UPN034' :  ['A*02:01',	'A*29:02',	'B*44:02',	'C*05:01',	'C*16:01'],
    'UPN035' :  ['A*01:01',	'A*26:01',	'B*07:02',	'B*40:01',	'C*02:02',	'C*07:02'],
    'UPN036' :  ['A*02:01',	'A*11:01',	'B*35:01',	'B*40:01',	'C*03:04',	'C*04:01'],
    'UPN037' :  ['A*02:01',	'A*24:02',	'B*15:01',	'B*44:02',	'C*03:04',	'C*16:01'],
    'UPN038' :  ['A*24:02',	'A*25:01',	'B*18:01',	'B*49:01',	'C*07:02',	'C*12:03'],
    'UPN039' :  ['A*03:01',	'A*24:02',	'B*35:01',	'C*04:01'],
    'UPN040' :  ['A*02:01',	'A*24:02',	'B*07:02',	'B*44:02',	'C*05:01',	'C*07:02'],
    'UPN041' :  ['A*02:01',	'A*11:01',	'B*07:02',	'B*15:01',	'C*03:04',	'C*07:02'],
    'UPN042' :  ['A*02:01',	'A*24:02',	'B*07:02',	'B*13:02',	'C*06:02',	'C*07:02'],
    'UPN043' :  ['A*01:01',	'A*68:01',	'B*44:02',	'B*55:01',	'C*03:04',	'C*05:01'],
    'UPN044' :  ['A*01:01',	'A*24:02',	'B*08:01',	'B*35:01',	'C*04:01',	'C*07:02'],
    'UPN045' :  ['A*03:01',	'A*24:02',	'B*51:01',	'C*14:02',	'C*16:01'],
    'UPN046' :  ['A*02:01',	'A*24:02',	'B*13:02',	'B*44:02',	'C*05:01',	'C*06:02'],
    'UPN047' :  ['A*02:01',	'A*26:01',	'B*07:02',	'B*40:01',	'C*03:04',	'C*07:02'],
    'UPN048' :  ['A*02:01',	'A*11:01',	'B*40:01',	'B*44:02',	'C*03:04',	'C*05:01'],
    'UPN049' :  ['A*02:01',	'B*08:01',	'B*40:01',	'C*03:04',	'C*07:02'],
    'UPN050' :  ['A*11:01',	'A*68:01',	'B*08:01',	'B*35:01',	'C*04:01',	'C*07:02'],
    'UPN051' :  ['A*01:01',	'A*24:02',	'B*15:01',	'B*40:01',	'C*02:02',	'C*07:02'],
    'UPN052' :  ['A*01:01',	'A*24:02',	'B*07:02',	'B*08:01',	'C*07:01',	'C*07:02'],
    'UPN053' :  ['A*02:01',	'A*30:01',	'B*13:02',	'B*35:01'],
    'UPN054' :  ['A*02:01',	'B*40:01',	'B*51:01',	'C*03:04',	'C*15:02'],
    'UPN055' :  ['A*02:01',	'B*55:01',	'B*57:01'],
    'UPN056' :  ['A*30:01',	'A*33:01',	'B*14:02',	'B*40:01',  'C*03:04',	'C*08:02'],
    'UPN057' :  ['A*02:01',	'A*03:01',	'B*27:05',	'B*57:01'],
    'UPN058' :  ['A*02:01',	'A*24:02',	'B*07:02',	'B*51:01',	'C*01:02',	'C*07:02'],
    'UPN059' :  ['A*03:01',	'A*26:01',	'B*07:02',	'B*38:01',	'C*07:02',	'C*12:03'],
    'UPN060' :  ['A*01:01',	'A*02:01',	'B*51:01',	'B*57:01',	'C*06:02',	'C*15:02'],
    'UPN061' :  ['A*02:01',	'A*11:01',	'B*39:06',	'C*06:02']
}

allele_map_II = { }

allele_map = allele_map_II

pin_file_dir = '/mnt/d/workspace/mhc-validator-2/data/PXD024871/XL/HLA-II'
mzml_folder = Path('/mnt/d/workspace/mhc-validator-2/data/PXD024871/XL/HLA-II')
output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/PXD024871/XL/HLA-II/all')
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
    if len(alleles) == 0 or 'CLL' in file_name:
        use_mhc_scores = False
        print('Warning: no allele specified for {}'.format(file_name))

    print(file_name)

    validator = MhcValidator(max_threads=20)  # Open a MHCvalidator instance, a new one has to be opened for each .pin file
    validator.set_mhc_params(alleles=alleles, mhc_class='II', min_pep_len=8)  # Load the alleles you specified above
    validator.load_data(pin, filetype='pin')  # Load the pin file
    validator.run(sequence_encoding=True, mhcflurry=use_mhc_scores, bigmhc=use_mhc_scores,
                  netmhcpan=use_mhc_scores, mixmhc2pred=use_mhc_scores,
                  autort=True, deeplc=False,
                  im2deep=False,
                  peptdeep=False,
                  # koina_predictors=rt_all_models + ms2_all_models + ccs_top_models,
                  koina_predictors=['Deeplc_hela_hf', 'ms2pip_Immuno_HCD', 'ms2pip_HCD2021'],
                  fine_tune=False,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  report_directory=output_folder / f'{file_name}')
