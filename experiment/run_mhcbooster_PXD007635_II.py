
import sys
from pathlib import Path

from mhcnames import normalize_allele_name

sys.path.append('/home/rw762/workspace/mhc-booster')
from mhcvalidator.validator import MhcValidator

allele_map_I = {
    'OvCa9_' :   ['A*02:01', 'A*03:01', 'B*07:02', 'B*40:02', 'C*07:02', 'C*12:01'],
    'OvCa10_' :  ['A*02:01', 'A*11:01', 'B*44:05', 'B*51:01', 'C*02:02', 'C*15:02'],
    'OvCa12_' :  ['A*24:02', 'A*31:01', 'B*35:03', 'B*49:01', 'C*07:01', 'C*12:03'],
    'OvCa13_' :  ['A*02', 'B*35', 'B*40', 'C*03', 'C*04'],   ###
    'OvCa15_' :  ['A*11:01', 'A*24:02', 'B*07:02', 'B*55:01', 'C*03:03', 'C*07:02'],
    'OvCa16_' :  ['A*02', 'B*40', 'B*4402', 'C*03', 'C*05'],   ###
    'OvCa23_' :  ['A*01', 'A*03', 'B*08', 'B*35', 'C*04', 'C*07'],   ###
    'OvCa28_' :  ['A*01:01', 'A*02:01', 'B*27:05', 'B*52:01', 'C*01:02', 'C*02:02'],
    'OvCa39_' :  ['A*25:01', 'A*31:01', 'B*07:02', 'B*18:01', 'C*12:03', 'C*07:02'],
    'OvCa41_' :  ['A*02', 'A*2402', 'B*18', 'B*51', 'C*02', 'C*12'],   ###
    'OvCa43_' :  ['A*02', 'A*32', 'B*18', 'B*35', 'C*04', 'C*07'],   ###
    'OvCa45_' :  ['A*01', 'A*23', 'B*08', 'B*44', 'C*04', 'C*07'],   ###
    'OvCa48_' :  ['A*02:01', 'A*25:01', 'B*15:01', 'B*41:02', 'C*03:04', 'C*17:01'],
    'OvCa53_' :  ['A*02', 'A*03', 'B*27', 'B*35', 'C*02', 'C*04'],   ###
    'OvCa54_' :  ['A*0201', 'A*11:01', 'B*35:01', 'B*35:03', 'C*04:01', 'C*12:03'],
    'OvCa57_' :  ['A*25', 'A*32', 'B*15', 'B*18', 'C*03', 'C*12'],   ###
    'OvCa58_' :  ['A*02', 'A*03', 'B*35', 'C*03', 'C*04'],   ###
    'OvCa59_' :  ['A*03', 'A*30', 'B*13', 'C*06'],   ###
    'OvCa60_' :  ['A*24:02', 'A*25:01', 'B*13:02', 'B*18:01', 'C*12:03', 'C*06:02'],
    'OvCa64_' :  ['A*01', 'A*25', 'B*08', 'C*07'],   ###
    'OvCa65_' :  ['A*01', 'A*2402', 'B*15', 'B*35', 'C*04', 'C*14'],   ###
    'OvCa66_' :  ['A*11:01', 'A*29:02', 'B*18:01', 'B*44:03', 'C*05:01', 'C*16:01'],
    'OvCa68_' :  ['A*02:01', 'A*01:01', 'B*44:02', 'B*37:01', 'C*06:02', 'C*05:01'],
    'OvCa69_' :  [],
    'OvCa70_' :  ['A*01', 'A*02', 'B*07', 'C*07'],   ###
    'OvCa72_' :  ['A*03:01', 'A*01:01', 'B*08:01', 'B*07:02', 'C*07:02', 'C*07:01'],
    'OvCa73_' :  ['A*01:01', 'B*08:01', 'C*07:01' ],
    'OvCa74_' :  ['A*02:01', 'B*18:01', 'B*51:01', 'C*07:02', 'C*15:02'],
    'OvCa79_' :  ['A*01:01', 'A*31:01', 'B*08:01', 'B*51:01', 'C*07:01', 'C*15:02'],
    'OvCa80_' :  ['A*25:01', 'A*32:01', 'B*18:01', 'B*39:01', 'C*12:03'],
    'OvCa81_' :  ['A*02:01', 'B*45:01', 'B*56:01', 'C*07:02', 'C*01:02'],
    'OvCa82_' :  ['A*01:01', 'A*03:01', 'B*08:01', 'B*38:01', 'C*07:01', 'C*12:03'],
    'OvCa83_' :  ['A*02', 'A*11', 'B*51', 'B*55', 'C*03', 'C*15'],   ###
    'OvCa84_' :  ['A*02:01', 'B*07:02', 'B*44:02', 'C*07:02', 'C*05:01'],
    'OvCa99_' :  ['A*02:01', 'A*24:02', 'B*13:02', 'B*40:01', 'C*03:04', 'C*06:02'],
    'OvCa100_' : ['A*02:01', 'B*07:02', 'B*41:02', 'C*07:02', 'C*17:01'],
    'OvCa103_' : ['A*02:01', 'A*24:02', 'B*27:02', 'B*27:05', 'C*02:02'],
    'OvCa104_' : ['A*03:01', 'B*07:02', 'B*35:08', 'C*04:01', 'C*07:02'],
    'OvCa105_' : ['A*26:01', 'A*68:01', 'B*18:01', 'B*55:01', 'C*03:03',' C*07:01'],
    'OvCa107_' : ['A*02:05', 'A*25:01', 'B*07:02', 'B*14:02', 'C*07:02',' C*08:02'],
    'OvCa109_' : ['A*02:01', 'A*23:01', 'B*40:01', 'B*49:01', 'C*07:01',' C*03:04'],
    'OvCa111_' : ['A*01:01', 'A*25:01', 'B*08:01', 'B*44:02', 'C*05:01',' C*07:01'],
    'OvCa114_' : ['A*29:02', 'B*44:03', 'C*16:01']
}

allele_map_II = {
    'OvCa9_' :   ['DQB1*03:01', 'DQA1*03:01', 'DQA1*05:01', 'DRB1*11:01', 'DRB1*04:01', 'DRB3*02:02', 'DRB4*01:01', 'DPB1*02:01', 'DPB1*13:01'],
    'OvCa10_' :  ['DQB1*02:02', 'DQB1*05:01', 'DQA1*01:01', 'DQA1*03:01', 'DRB1*01:01', 'DRB1*09:01', 'DRB4*01:01', 'DPB1*04:01', 'DPB1*05:01'],
    'OvCa12_' :  ['DQB1*03:01', 'DQB1*05:04', 'DQA1*01:02', 'DQA1*03:01', 'DRB1*01:01', 'DRB1*04:01', 'DRB4*01:01', 'DPB1*02:01', 'DPB1*05:01'],
    'OvCa13_' :  ['DQB1*04', 'DQB1*06', 'DRB1*08', 'DRB1*13'],
    'OvCa15_' :  ['DQB1*03:01', 'DQA1*05:01', 'DRB1*11:01', 'DRB1*03:17', 'DRB3*02:02', 'DPB1*03:01'],
    'OvCa16_' :  ['DQB1*06', 'DRB1*08', 'DRB1*13', 'DRB1*14', 'DRB3'],
    'OvCa23_' :  ['DQB1*02', 'DQB1*03', 'DRB1*03', 'DRB1*12', 'DRB3'],
    'OvCa28_' :  ['DQB1*05:01', 'DQB1*06:01', 'DQA1*01:01', 'DQA1*03:01', 'DRB1*01:03', 'DRB1*15:02', 'DRB5*01:02', 'DPB1*04:01'],
    'OvCa39_' :  ['DQB1*06:02', 'DQA1*01:02', 'DRB1*15:01', 'DRB1*16:09', 'DRB5*01:01', 'DRB5*01:11', 'DPB1*04:01', 'DPB1*04:02'],
    'OvCa41_' :  ['DQB1*03', 'DQ7', 'DRB1*11', 'DRB3'],
    'OvCa43_' :  ['DQB1*03', 'DQB1*05', 'DQ9', 'DRB1*01', 'DRB1*07', 'DRB4'],
    'OvCa45_' :  ['DQB1*02', 'DRB1*03', 'DRB1*07', 'DRB3', 'DRB4'],
    'OvCa48_' :  ['DQB1*03:02', 'DQB1*03:04', 'DQA1*03:01', 'DRB1*04:01', 'DRB1*13:03', 'DRB3*01:01', 'DRB4*01:01', 'DPB1*02:01'],
    'OvCa53_' :  ['DQB1*02', 'DQB1*03', 'DQ7', 'DRB1*03', 'DRB1*11', 'DRB3'],
    'OvCa54_' :  ['DQB1*05:01', 'DQB1*05:03', 'DQA1*01:01', 'DRB1*01:03', 'DRB1*14:01', 'DRB3*02:02', 'DPB1*04:01', 'DPB1*02:01'],
    'OvCa57_' :  ['DQB1*05', 'DQB1*06', 'DRB1*01', 'DRB1*15', 'DRB5'],
    'OvCa58_' :  ['DQB1*05', 'DRB1*01'],
    'OvCa59_' :  ['DQB1*02', 'DRB1*07', 'DRB4'],
    'OvCa60_' :  ['DRB1*08:01', 'DRB1*13:01', 'DQB1*04:02', 'DQB1*06:03', 'DQA1*04:01', 'DQA1*01:03', 'DPB1*02:01', 'DPB1*03:01'],
    'OvCa64_' :  ['DQB1*02', 'DRB1*03', 'DRB3'],
    'OvCa65_' :  ['DQB1*03', 'DQB1*05', 'DRB1*10', 'DRB1*11', 'DRB3'],
    'OvCa66_' :  ['DRB1*03', 'DRB*07:01', 'DRB3*02:02', 'DRB4*01:01', 'DQB1*02:01', 'DQB1*02:02', 'DQA1*02:01', 'DQA1*05:01', 'DPB1*02:02', 'DPB1*03:01'],
    'OvCa68_' :  ['DRB1*10:01', 'DRB1*04:01', 'DRB4*04:01', 'DQB1*05:01', 'DQB1*03:01', 'DQA1*01:01', 'DPB1*04:01'],
    'OvCa69_' :  [],
    'OvCa70_' :  ['DQB1*03', 'DQB1*05', 'DRB1*09', 'DRB1*14', 'DRB3', 'DRB4'],
    'OvCa72_' :  ['DRB1*01:01', 'DRB1*03:01', 'DRB3*01:01', 'DQB1*05:01', 'DQB1*02:01', 'DQA1*01:01', 'DPB1*04:01'],
    'OvCa73_' :  ['DRB1*03:01', 'DRB1*03:42', 'DRB3*01:01', 'DRB3*01:14', 'DQB1*02:01', 'DQA1*05:01', 'DPB1*04:01'],
    'OvCa74_' :  ['DRB1*11:04', 'DRB1*07:01', 'DRB3*02:02', 'DRB4*01:01', 'DQB1*03:01', 'DQB1*02:02', 'DQA1*02:01', 'DQA1*05:01', 'DPB1*04:02', 'DPB1*02:01'],
    'OvCa79_' :  ['DQB1*03:03', 'DQA1*02:01', 'DRB1*07:01', 'DRB1*09:01', 'DRB4*01:01', 'DPB1*13:01', 'DPB1*02:01'],
    'OvCa80_' :  ['DRB1*01:01', 'DRB1*12:01', 'DRB3*02:02', 'DQB1*03:01', 'DQB1*05:01', 'DQA1*01:01', 'DQA1*05:01', 'DPB1*04:01'],
    'OvCa81_' :  ['DRB1*04:02', 'DRB1*11:01', 'DRB4*01:01', 'DRB3*02:02', 'DQB1*03:01', 'DQB1*03:02'],
    'OvCa82_' :  ['DRB1*04:02', 'DRB1*03:01', 'DRB4*01:01', 'DRB3*01:01', 'DQB1*02:01', 'DQB1*03:02', 'DQA1*03:01', 'DQA1*05:01', 'DPB1*04:01', 'DPB1*13:01'],
    'OvCa83_' :  ['DQB1*03', 'DQB1*05', 'DRB1*09', 'DRB1*14', 'DRB3', 'DRB4'],
    'OvCa84_' :  ['DRB1*15:01', 'DRB5*01:01', 'DQB1*06:02', 'DQA1*01:02', 'DPB1*04:01', 'DPB1*04:02'],
    'OvCa99_' :  ['DRB1*04:01', 'DRB1*07:01', 'DRB4*01:01', 'DQB1*02:02', 'DQB1*03:01', 'DQA1*02:01', 'DQA1*03:01', 'DPB1*02:01', 'DPB1*04:01'],
    'OvCa100_' : ['DRB1*13:03', 'DRB1*15:01', 'DRB3*01:01', 'DRB5*01:01', 'DQB1*03:01', 'DQB1*06:02', 'DQA1*01:02', 'DQA1*05:01', 'DPB1*03:01', 'DPB1*05:01'],
    'OvCa103_' : ['DRB1*01:01', 'DRB1*15:01', 'DRB5*01:01', 'DQB1*05:01', 'DQB1*05:02', 'DQA1*01:01', 'DQA1*01:02', 'DPB1*02:01', 'DPB1*04:02'],
    'OvCa104_' : ['DRB1*03:01', 'DRB3*02:02', 'DQB1*02:01', 'DQB1*03:01', 'DQA1*05:01', 'DPB1*04:02'],
    'OvCa105_' : ['DRB1*13:02', 'DRB1*15:01', 'DRB3*03:01', 'DRB5*01:01', 'DQB1*06:03', 'DQB1*06:04', 'DQA1*01:02', 'DPB1*04:01'],
    'OvCa107_' : ['DQB1*03:01', 'DQ7', 'DRB1*01:02', 'DRB1*11:06', 'DRB3'],
    'OvCa109_' : ['DRB1*01:01', 'DRB1*13:02', 'DQB1*05:01', 'DQB1*06:04'],
    'OvCa111_' : ['DRB1*01:01', 'DRB1*03:01', 'DQB1*02:01', 'DQB1*05:01'],
    'OvCa114_' : ['DRB1*07:01', 'DQB1*02:02'],
}

allele_map = allele_map_II

pin_file_dir = '/mnt/d/workspace/mhc-validator-2/data/PXD007635/HLA-II'
mzml_folder = Path('/mnt/d/workspace/mhc-validator-2/data/PXD007635/HLA-II')
output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/PXD007635/HLA-II/test')
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
    if len(alleles) == 0:
        use_mhc_scores = False
        print('Warning: no allele specified for {}'.format(file_name))

    print(file_name)

    validator = MhcValidator(max_threads=20)  # Open a MHCvalidator instance, a new one has to be opened for each .pin file
    validator.set_mhc_params(alleles=alleles, mhc_class='II', min_pep_len=9)  # Load the alleles you specified above
    validator.load_data(pin, filetype='pin')  # Load the pin file
    validator.run(sequence_encoding=True, mhcflurry=use_mhc_scores, bigmhc=use_mhc_scores,
                  netmhcpan=use_mhc_scores, mixmhc2pred=use_mhc_scores,
                  autort=False, deeplc=False,
                  im2deep=False,
                  peptdeep=False,
                  # koina_predictors=rt_all_models + ms2_all_models + ccs_top_models,
                  koina_predictors=['Prosit_2024_irt_cit', 'Prosit_2019_irt', 'ms2pip_Immuno_HCD', 'ms2pip_HCD2021'],
                  fine_tune=False,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  report_directory=output_folder / f'{file_name}')
