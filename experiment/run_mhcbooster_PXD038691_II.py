
import sys
from pathlib import Path

sys.path.append('/home/rw762/workspace/mhc-booster')
from mhcvalidator.validator import MhcValidator

allele_map_I = {
    'UPN01' :   ['A*11:01', 'B*07:02', 'B*44:03', 'C*07:02', 'C*16:01'],
    'UPN02' :   ['A*01:01', 'B*07:02', 'B*08:01', 'C*07:01', 'C*07:02'],
    'UPN03' :   ['A*03:01', 'A*23:01', 'B*18:01', 'B*44:02', 'C*05:01', 'C*12:03'],
    'UPN04' :   ['A*02:01', 'A*31:01', 'B*40:01', 'B*44:02', 'C*03:04', 'C*05:01'],
    'UPN05' :   ['A*01:01', 'A*29:02', 'B*08:01', 'B*44:03', 'C*07:01', 'C*16:01'],
    'UPN06' :   ['A*03:01', 'A*23:01', 'B*08:01', 'B*15:01', 'C*03:03', 'C*07:01'],
    'UPN07' :   ['A*02:01', 'A*24:10', 'B*14:02', 'B*35:30', 'C*04:01', 'C*08:02'],
    'UPN08' :   ['A*01:01', 'A*25:01', 'B*18:01', 'B*57:01', 'C*06:02', 'C*12:03'],
    'UPN09' :   ['A*01:01', 'B*08:01', 'B*13:02', 'C*06:02', 'C*07:01'],
    'UPN10' :   ['A*02:01', 'A*23:01', 'B*44:02', 'B*44:03', 'C*02:02', 'C*04:01'],
    'UPN11' :   ['A*01:01', 'A*03:01', 'B*08:01', 'B*14:02', 'C*07:01', 'C*08:02'],
    'UPN12' :   ['A*02:01', 'A*02:05', 'B*51:01', 'B*58:01', 'C*03:02', 'C*14:02'],
    'UPN13' :   ['A*11:01', 'A*68:01', 'B*15:01', 'B*51:01', 'C*03:03', 'C*14:02'],
    'UPN14' :   ['A*02:01', 'A*11:01', 'B*35:01', 'B*44:02', 'C*04:01', 'C*05:01'],
    'UPN15' :   ['A*02:01', 'A*11:01', 'B*08:01', 'B*57:01', 'C*07:01'],
    'UPN16' :   ['A*03:01', 'A*26:01', 'B*07:02', 'B*35:01', 'C*04:01', 'C*07:02'],
    'UPN17' :   ['A*03:01', 'A*66:01', 'B*14:01', 'B*41:02', 'C*08:02', 'C*17:01'],
    'UPN18' :   ['A*02:01', 'A*24:02', 'B*07:02', 'B*51:01', 'C*07:02', 'C*14:02'],
    'UPN19' :   ['A*02:01', 'A*31:01', 'B*07:02', 'C*07:02'],
    'UPN20' :   ['A*02:05', 'A*24:02', 'B*07:02', 'B*55:01', 'C*01:02', 'C*07:02'],
    'UPN21' :   ['A*02:01', 'A*25:01', 'B*44:02', 'B*57:01', 'C*05:01', 'C*06:02'],
    'UPN22' :   ['A*01:01', 'A*03:01', 'B*07:02', 'B*08:01', 'C*07:01', 'C*07:02'],
    'UPN23' :   ['A*02:02', 'A*11:01', 'B*35:02', 'B*44:05', 'C*02:02', 'C*04:01'],
    'UPN24' :   ['A*03:01', 'A*24:02', 'B*07:02', 'B*51:01', 'C*07:02'],
    'UPN25' :   ['A*26:01', 'A*32:01', 'B*15:01', 'B*38:01', 'C*03:03', 'C*12:03'],
    'UPN26' :   ['A*02:01', 'B*35:03', 'B*40:01', 'C*03:04', 'C*04:01'],
    'UPN27' :   ['A*02:01', 'A*33:01', 'B*13:02', 'B*14:02', 'C*06:02', 'C*08:02'],
    'UPN28' :   ['A*30:01', 'A*33:03', 'B*13:02', 'B*44:02', 'C*06:02', 'C*07:04'],
    'UPN29' :   ['A*02:01', 'A*68:02', 'B*15:01', 'B*53:01', 'C*03:03', 'C*04:01'],
    'UPN30' :   ['A*24:02', 'A*25:01', 'B*07:02', 'B*52:01', 'C*07:02', 'C*12:02'],
    'UPN31' :   ['A*24:02', 'A*68:01', 'B*51:01', 'B*56:01', 'C*01:02', 'C*02:02'],
    'UPN32' :   ['A*01:01', 'A*24:02', 'B*35:01', 'B*57:01', 'C*04:01', 'C*06:02'],
    'UPN33' :   ['A*02:01', 'A*26:01', 'B*07:10', 'B*40:01', 'C*03:04', 'C*07:02'],
    'UPN34' :   ['A*02:01', 'B*40:01', 'C*03:04'],
    'UPN35' :   ['A*01:01', 'A*31:01', 'B*08:01', 'B*55:01', 'C*03:03', 'C*07:01'],
    'UPN36' :   ['A*02:01', 'A*11:01', 'B*14:01', 'B*55:01', 'C*03:03', 'C*08:02'],
    'UPN37' :   ['A*01:01', 'A*02:01', 'B*08:01', 'B*15:01', 'C*03:04', 'C*07:01'],
    'UPN38' :   ['A*02:01', 'A*31:01', 'B*27:05', 'B*44:29', 'C*02:02', 'C*07:02'],
    'UPN39' :   ['A*01:01', 'A*03:01', 'B*07:02', 'B*35:03', 'C*06:02', 'C*07:02'],
    'UPN40' :   ['A*02:01', 'A*23:01', 'B*44:02', 'B*49:01', 'C*05:01', 'C*07:01'],
    'UPN41' :   ['A*11:01', 'A*24:02', 'B*49:01', 'B*52:01', 'C*07:01', 'C*12:02'],
    'UPN42' :   ['A*01:01', 'A*24:02', 'B*07:02', 'B*08:01', 'C*07:01', 'C*07:02'],
    'UPN43' :   ['A*01:01', 'A*68:01', 'B*35:03', 'B*52:01', 'C*04:01', 'C*12:02'],
    'UPN44' :   ['A*01:01', 'A*02:01', 'B*08:01', 'B*38:01', 'C*07:01', 'C*12:03'],
    'UPN45' :   ['A*01:01', 'A*29:02', 'B*44:03', 'B*57:01', 'C*06:02', 'C*16:01'],
    'UPN46' :   ['A*02:01', 'A*31:01', 'B*15:01', 'B*18:01', 'C*03:04', 'C*07:01'],
    'UPN47' :   ['A*01:01', 'A*02:01', 'B*08:01', 'B*44:02', 'C*05:01', 'C*07:01'],
    'UPN48' :   ['A*24:02', 'A*68:01', 'B*15:17', 'B*40:01', 'C*03:04', 'C*07:01'],
    'UPN49' :   ['A*02:01', 'A*03:01', 'B*27:02', 'B*35:03', 'C*02:02', 'C*04:01'],
    'UPN50' :   ['A*03:01', 'A*30:01', 'B*35:01', 'B*44:03', 'C*04:01', 'C*16:01'],
    'UPN51' :   ['A*03:01', 'A*32:03', 'B*15:17', 'B*35:01', 'C*04:01', 'C*15:02'],
    'UPN52' :   ['A*01:01', 'A*02:06', 'B*27:05', 'B*35:01', 'C*03:03', 'C*04:01']

}

allele_map_II = {
    'UPN01':    ['DRB1*03:01', 'DRB1*07:01', 'DQB1*02:01', 'DQB1*02:02', 'DQA1*02:01', 'DQA1*05:01'],
    'UPN02':    ['DRB1*07:01', 'DRB1*15:01', 'DQB1*02:02', 'DQB1*06:02', 'DQA1*01:02', 'DQA1*02:01'],
    'UPN03':    ['DRB1*13:01', 'DRB1*14:01', 'DQB1*05:03', 'DQB1*06:03', 'DQA1*01:01', 'DQA1*01:03'],
    'UPN04':    ['DRB1*04:01', 'DRB1*09:01', 'DQB1*03:02', 'DQB1*03:03', 'DQA1*03:01'],
    'UPN05':    ['DRB1*01:01', 'DRB1*07:01', 'DQB1*02:02', 'DQB1*05:01', 'DQA1*01:02', 'DQA1*02:01'],
    'UPN06':    ['DRB1*11:01', 'DQB1*06:02', 'DQA1*01:01', 'DQA1*05:05'],
    'UPN07':    ['DRB1*08:01', 'DRB1*12:04', 'DQB1*04:02', 'DQB1*05:02', 'DQA1*01:01', 'DQA1*04:01'],
    'UPN08':    ['DRB1*07:01', 'DRB1*08:01', 'DQB1*03:03', 'DQB1*04:02', 'DQA1*02:01', 'DQA1*04:01'],
    'UPN09':    ['DRB1*03:01', 'DRB1*07:01', 'DQB1*02:01', 'DQB1*02:02', 'DQA1*02:01', 'DQA1*05:01'],
    'UPN10':    ['DRB1*07:01', 'DQB1*02:02', 'DQA1*02:02'],
    'UPN11':    ['DRB1*03:01', 'DRB1*15:01', 'DQB1*02:01', 'DQB1*06:02', 'DQA1*01:01', 'DQA1*05:01'],
    'UPN12':    ['DRB1*03:01', 'DRB1*16:01', 'DQB1*02:01', 'DQB1*05:02', 'DQA1*01:01', 'DQA1*05:01'],
    'UPN13':    ['DRB1*13:01', 'DRB1*15:01', 'DQB1*06:02', 'DQB1*06:03', 'DQA1*01:02', 'DQA1*01:03'],
    'UPN14':    ['DRB1*03:01', 'DRB1*04:01', 'DQB1*03:02', 'DQB1*05:01', 'DQA1*01:01', 'DQA1*03:01'],
    'UPN15':    ['DRB1*03:01', 'DRB1*07:01', 'DQB1*02:01', 'DQB1*03:03', 'DQA1*02:01', 'DQA1*05:01'],
    'UPN16':    ['DRB1*01:01', 'DRB1*04:04', 'DQB1*03:02', 'DQB1*05:01'],
    'UPN17':    ['DRB1*13:02', 'DRB1*13:03', 'DQB1*03:01', 'DQB1*06:09', 'DQA1*01:02', 'DQA1*05:05'],
    'UPN18':    ['DRB1*11:04', 'DRB1*13:01', 'DQB1*03:01', 'DQB1*05:01', 'DQA1*01:03', 'DQA1*05:05'],
    'UPN19':    ['DRB1*14:01', 'DRB1*15:01', 'DQB1*05:03', 'DQB1*06:02', 'DQA1*01:01', 'DQA1*01:02'],
    'UPN20':    ['DRB1*13:01', 'DRB1*14:01', 'DQB1*05:03', 'DQB1*06:03', 'DQA1*01:01', 'DQA1*01:03'],
    'UPN21':    ['DRB1*13:01', 'DRB1*13:03', 'DQB1*03:01', 'DQB1*06:03', 'DQA1*01:01', 'DQA1*05:03'],
    'UPN22':    ['DRB1*03:01', 'DRB1*13:02', 'DQB1*02:01', 'DQB1*06:04', 'DQA1*01:02', 'DQA1*05:01'],
    'UPN23':    ['DRB1*03:01', 'DRB1*16:01', 'DQB1*02:01', 'DQB1*05:02', 'DQA1*02:02', 'DQA1*05:01'],
    'UPN24':    ['DRB1*04:07', 'DRB1*15:01', 'DQB1*03:01', 'DQB1*06:02', 'DQA1*01:02', 'DQA1*02:01'],
    'UPN25':    ['DRB1*04:04', 'DRB1*12:01', 'DQB1*03:01', 'DQB1*03:02', 'DQA1*03:01', 'DQA1*05:05'],
    'UPN26':    ['DRB1*13:01', 'DRB1*14:54', 'DQB1*05:03', 'DQB1*06:03'],
    'UPN27':    ['DRB1*01:02', 'DRB1*07:01', 'DQB1*02:02', 'DQB1*05:01', 'DPB1*04:01'],
    'UPN28':    ['DRB1*07:01', 'DRB1*11:01', 'DQB1*02:02', 'DQB1*03:01', 'DQA1*02:01', 'DQA1*05:05'],
    'UPN29':    ['DRB1*11:04', 'DRB1*13:02', 'DQB1*03:01', 'DQB1*06:04', 'DPB1*03:01', 'DPB1*04:02'],
    'UPN30':    ['DRB1*07:01', 'DRB1*15:02', 'DQB1*03:03', 'DQB1*06:01'],
    'UPN31':    ['DRB1*04:01', 'DRB1*08:01', 'DQB1*03:01'],
    'UPN32':    ['DRB1*01:01', 'DRB1*07:01', 'DQB1*03:03', 'DQB1*05:01'],
    'UPN33':    ['DRB1*03:01', 'DRB1*13:02', 'DQB1*02:01', 'DQB1*06:04', 'DPB1*04:01'],
    'UPN34':    ['DRB1*08:01', 'DRB1*13:02', 'DQB1*04:02', 'DQB1*06:04', 'DPB1*03:01', 'DPB1*04:02'],
    'UPN35':    ['DRB1*14:01', 'DRB1*14:06', 'DQB1*03:01', 'DQB1*05:03', 'DQA1*01:01', 'DQA1*05:03'],
    'UPN36':    ['DRB1*14:54', 'DRB1*15:01', 'DQB1*05:03', 'DQB1*06:02', 'DPB1*04:02', 'DPB1*06:01'],
    'UPN37':    ['DRB1*04:01', 'DQB1*03:02', 'DPB1*04:01'],
    'UPN38':    ['DRB1*13:01', 'DQB1*06:03', 'DQA1*01:03'],
    'UPN39':    ['DRB1*13:05', 'DRB1*15:01', 'DQB1*03:01', 'DQB1*06:02', 'DQA1*01:01', 'DQA1*05:05'],
    'UPN40':    ['DRB1*04:01', 'DRB1*15:01', 'DQB1*03:02', 'DQB1*06:02', 'DQA1*01:02', 'DQA1*03:01'],
    'UPN41':    ['DRB1*01:01', 'DRB1*11:04', 'DQB1*03:01', 'DQB1*05:01', 'DQA1*01:01', 'DQA1*05:05'],
    'UPN42':    [],
    'UPN43':    ['DRB1*11:03', 'DRB1*15:02', 'DQB1*03:01', 'DQB1*06:01'],
    'UPN44':    ['DRB1*03:01', 'DRB1*13:01', 'DQB1*02:01', 'DQB1*06:03'],
    'UPN45':    ['DRB1*07:01', 'DRB1*14:01', 'DQB1*02:02', 'DQB1*05:03'],
    'UPN46':    ['DRB1*04:01', 'DRB1*13:01', 'DQB1*03:02', 'DQB1*06:03', 'DPB1*04:01'],
    'UPN47':    ['DRB1*03:01', 'DRB1*14:54', 'DQB1*02:01', 'DQB1*05:03', 'DPB1*01:01', 'DPB1*09:01'],
    'UPN48':    ['DRB1*13:02', 'DQB1*06:04', 'DQA1*01:01', 'DQA1*01:03'],
    'UPN49':    ['DRB1*01:01', 'DRB1*12:01', 'DQB1*03:01', 'DQB1*05:01', 'DPB1*04:01', 'DPB1*13:01'],
    'UPN50':    ['DRB1*07:01', 'DRB1*11:01', 'DQB1*02:02', 'DQB1*03:01'],
    'UPN51':    ['DRB1*01:01', 'DQB1*05:01'],
    'UPN52':    ['DRB1*15:01', 'DQB1*06:02']
}

allele_map = allele_map_II

pin_file_dir = '/mnt/d/workspace/mhc-validator-2/data/PXD038691/HLA-II'
mzml_folder = Path('/mnt/d/workspace/mhc-validator-2/data/PXD038691/HLA-II')
output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/PXD038691/HLA-II/test')
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
                  # koina_predictors=ms2_all_models + rt_all_models,
                  koina_predictors=['Prosit_2024_irt_cit', 'Deeplc_hela_hf', 'AlphaPeptDeep_rt_generic', 'ms2pip_HCD2021', 'ms2pip_Immuno_HCD'],
                  fine_tune=False,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  report_directory=output_folder / f'{file_name}')
