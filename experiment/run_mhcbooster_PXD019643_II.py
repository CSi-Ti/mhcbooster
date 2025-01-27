
import sys
from pathlib import Path

sys.path.append('/home/rw762/workspace/mhc-booster')
from mhcvalidator.validator import MhcValidator

allele_map_I = {
    'AUT01-DN02':   ['HLA-A1101', 'HLA-A6801', 'HLA-B1501', 'HLA-B3503', 'HLA-C0303', 'HLA-C0401'], # 21
    'AUT01-DN03':   ['HLA-A0101', 'HLA-A1101', 'HLA-B1501', 'HLA-B3501', 'HLA-C0303', 'HLA-C0401'], # 36
    'AUT01-DN04':   ['HLA-A0201', 'HLA-A2301', 'HLA-B2705', 'HLA-B5001', 'HLA-C0202', 'HLA-C0602'], # 39
    'AUT01-DN05':   ['HLA-A0101', 'HLA-A1101', 'HLA-B0702', 'HLA-B4901', 'HLA-C0701', 'HLA-C0702'], # 21
    'AUT01-DN06':   ['HLA-A0301', 'HLA-A6802', 'HLA-B0702', 'HLA-B1402', 'HLA-C0702', 'HLA-C0802'], # 66
    'AUT01-DN08':   ['HLA-A3201', 'HLA-A6801', 'HLA-B1501', 'HLA-B4402', 'HLA-C0303', 'HLA-C0704'], # 66
    'AUT01-DN09':   ['HLA-A2402', 'HLA-A3001', 'HLA-B1302', 'HLA-B3508', 'HLA-C0401', 'HLA-C0602'], # 63
    'AUT01-DN11':   ['HLA-A0101', 'HLA-A6901', 'HLA-B3701', 'HLA-B4901', 'HLA-C0602', 'HLA-C0701'], # 75
    'AUT01-DN12':   ['HLA-A0201', 'HLA-A1101', 'HLA-B1501', 'HLA-B3501', 'HLA-C0304', 'HLA-C0401'], # 72
    'AUT01-DN13':   ['HLA-A0205', 'HLA-A1101', 'HLA-B4002', 'HLA-B5801', 'HLA-C0202', 'HLA-C0701'], # 72
    'AUT01-DN14':   ['HLA-A0201', 'HLA-A6802', 'HLA-B1402', 'HLA-B2705', 'HLA-C0202', 'HLA-C0802'], # 66
    'AUT01-DN15':   ['HLA-A0101', 'HLA-A0201', 'HLA-B0801', 'HLA-B4402', 'HLA-C0701', 'HLA-C0704'], # 51
    'AUT01-DN16':   ['HLA-A0101', 'HLA-A2402', 'HLA-B0801', 'HLA-B4101', 'HLA-C0701', 'HLA-C1701'], # 15
    'AUT01-DN17':   ['HLA-A0301', 'HLA-A2402', 'HLA-B3503', 'HLA-B4501', 'HLA-C0401', 'HLA-C1601'], # 66
    'OVA01-DN278':  ['HLA-A0201', 'HLA-B4402', 'HLA-C0501'], # 15
    'OVA01-DN281':  ['HLA-A1101', 'HLA-A2601', 'HLA-B0801', 'HLA-B3501', 'HLA-C0702', 'HLA-C0401'], # 15
    'THY01-DN1':    ['HLA-A0301', 'HLA-A2902', 'HLA-B0702', 'HLA-B4403', 'HLA-C0702', 'HLA-C1601'], # 3
    'THY01-DN3':    ['HLA-A2402', 'HLA-A2501', 'HLA-B1801', 'HLA-B4101', 'HLA-C1203', 'HLA-C1701'], # 3
    'THY01-DN4':    ['HLA-A0201', 'HLA-A2608', 'HLA-B1501', 'HLA-B4402', 'HLA-C0304', 'HLA-C0501'], # 3
    'THY01-DN5':    ['HLA-A0101', 'HLA-A0301', 'HLA-B0706', 'HLA-B0702', 'HLA-C0702', 'HLA-C1505'], # 3
    'THY01-DN6':    ['HLA-A0101', 'HLA-A2501', 'HLA-B1302', 'HLA-B3901', 'HLA-C0602', 'HLA-C1203'], # 3
}

allele_map_II = {
    'AUT01-DN02':   ['DRB1*04:01', 'DRB4*01:01', 'DQA1*03:01', 'DQB1*03:01', 'DPB1*04:01', 'DPB1*04:02'],
    'AUT01-DN03':   ['DRB1*07:01', 'DRB1*11:03', 'DRB3*02:02', 'DRB4*01:01', 'DQB1*02:02', 'DQB1*03:01', 'DQA1*02:01', 'DQA1*05:01', 'DPB1*03:01', 'DPB1*04:01'],
    'AUT01-DN04':   ['DRB1*13:01', 'DRB1*07:01', 'DRB3*01:01', 'DRB4*01:01', 'DQB1*02:01', 'DQB1*06:03', 'DQA1*01:03', 'DQA1*02:01', 'DPB1*03:01', 'DPB1*04:01'],
    'AUT01-DN05':   ['DRB1*11:01', 'DRB1*14:01', 'DRB3*02:02', 'DQB1*03:01', 'DQB1*05:03', 'DQA1*01:01', 'DQA1*05:01', 'DPB1*04:01'],
    'AUT01-DN06':   ['DRB1*13:03', 'DRB1*08:01', 'DRB3*01:01', 'DQB1*03:01', 'DQB1*04:02'],
    'AUT01-DN08':   ['DRB1*13:03', 'DRB1*14:01', 'DQB1*03:01', 'DQA1*05:05'],
    'AUT01-DN09':   ['DRB1*07:01', 'DQB1*02:02', 'DQA1*02:01'],
    'AUT01-DN11':   ['DRB1*07:01', 'DRB1*13:02', 'DQB1*03:03', 'DQB1*06:04', 'DQA1*01:02', 'DQA1*02:01'],
    'AUT01-DN12':   ['DRB1*01:01', 'DRB1*12:01', 'DQB1*03:01', 'DQB1*05:01', 'DQA1*01:01', 'DQA1*05:05'],
    'AUT01-DN13':   ['DRB1*10:01', 'DRB1*15:01', 'DQB1*05:01', 'DQB1*06:02', 'DQA1*01:01', 'DQA1*01:02'],
    'AUT01-DN14':   ['DRB1*04:01', 'DRB1*13:03', 'DQB1*03:01', 'DQB1*03:02', 'DQA1*03:02', 'DQA1*05:05'],
    'AUT01-DN15':   ['DRB1*03:01', 'DRB1*11:01', 'DQB1*02:01', 'DQB1*03:01', 'DQA1*05:01'],
    'AUT01-DN16':   ['DRB1*03:01', 'DRB1*04:05', 'DQB1*02:01', 'DQB1*02:02', 'DQA1*03:02', 'DQA1*05:01'],
    'AUT01-DN17':   ['DRB1*14:01', 'DRB1*15:01', 'DQB1*05:03', 'DQB1*06:02'],
    'OVA01-DN278':  ['DRB1*11:01', 'DRB1*13:01', 'DQB1*03:01', 'DQB1*06:03'],
    'OVA01-DN281':  ['DRB1*04:04', 'DRB1*11:08', 'DQB1*03:02', 'DQB1*03:01'],
    'THY01-DN1':    ['DRB1*01:01', 'DRB1*15:01', 'DRB5*01:01', 'DQB1*05:01', 'DQB1*06:02', 'DQA1*01:01', 'DQA1*01:02', 'DPB1*04:01', 'DPB1*14:01', 'DPA1*01:03', 'DPA1*02:01'],
    'THY01-DN3':    ['DRB1*04:05', 'DRB1*15:01', 'DRB4*01:03', 'DRB5*01:01', 'DQB1*02:02', 'DQB1*06:02', 'DQA1*01:02', 'DQA1*03:03', 'DPB1*02:01', 'DPB1*09:01', 'DPA1*01:03', 'DPA1*02:01'],
    'THY01-DN4':    ['DRB1*11:01', 'DRB1*15:01', 'DRB3*02:02', 'DRB5*01:01', 'DQB1*03:01', 'DQB1*06:02', 'DQA1*01:02', 'DQA1*05:05', 'DPB1*02:01', 'DPB1*04:01', 'DPA1*01:03'],
    'THY01-DN5':    ['DRB1*04:05', 'DRB1*15:01', 'DRB4*01:03', 'DRB5*01:01', 'DQB1*03:02', 'DQB1*06:02', 'DQA1*01:02', 'DQA1*03:03', 'DPB1*104:01', 'DPB1*04:01', 'DPA1*01:03'],
    'THY01-DN6':    ['DRB1*10:01', 'DRB1*16:01', 'DRB5*02:02', 'DQB1*05:01', 'DQB1*05:02', 'DQA1*01:05', 'DQA1*01:02', 'DPB1*02:01', 'DPB1*17:01', 'DPA1*01:03', 'DPA1*02:01'],
}
allele_map = allele_map_II
pin_file_dir = '/mnt/d/workspace/mhc-validator-2/data/PXD019643/HLA-II'
mzml_folder = Path('/mnt/d/workspace/mhc-validator-2/data/PXD019643/HLA-II')
output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/PXD019643/HLA-II/app')
output_folder.mkdir(parents=True, exist_ok=True)

pin_files = list(Path(pin_file_dir).rglob('*.pin'))

rt_all_models = ['Deeplc_hela_hf', 'AlphaPeptDeep_rt_generic', 'Prosit_2019_irt', 'Prosit_2024_irt_cit', 'Chronologer_RT']
ms2_all_models = ['Prosit_2019_intensity', 'Prosit_2024_intensity_cit', 'Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_CID', 'Prosit_2020_intensity_HCD',
                  'ms2pip_HCD2021', 'ms2pip_timsTOF2023', 'ms2pip_Immuno_HCD', 'ms2pip_timsTOF2024']

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
    validator.run(sequence_encoding=False, mhcflurry=use_mhc_scores, bigmhc=use_mhc_scores,
                  netmhcpan=use_mhc_scores, mixmhc2pred=use_mhc_scores,
                  autort=False, deeplc=False,
                  im2deep=False,
                  peptdeep=False,
                  # koina_predictors=rt_all_models + ms2_all_models,
                  # koina_predictors=['Prosit_2024_irt_cit', 'Prosit_2019_irt', 'Prosit_2020_intensity_CID', 'Prosit_2024_intensity_cit'],
                  fine_tune=False,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  report_directory=output_folder / f'{file_name}')

