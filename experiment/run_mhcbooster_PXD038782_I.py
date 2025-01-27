
import sys
from pathlib import Path

sys.path.append('/home/rw762/workspace/mhc-booster')
from mhcvalidator.validator import MhcValidator

allele_map_I = {
    'UDN01' :	['A*01', 'A*11', 'B*15', 'B*35', 'C*03', 'C*04'],
    'UDN02' :	['A*02', 'A*23', 'B*27', 'B*50', 'C*0202', 'C*0602'],
    'UDN03' :	['A*01', 'A*11', 'B*0702', 'B*49', 'C*07', 'C*0702'],
    'UDN04' :	['A*03', 'A*68', 'B*0702', 'B*14', 'C*07', 'C*08'],
    'UDN05' :	['A*2402', 'A*30', 'B*13', 'B*35', 'C*04', 'C*0602'],
    'UDN06' :	['A*01', 'A*69', 'B*37', 'B*49', 'C*0602', 'C*07'],
    'UDN07' :	['A*02', 'A*11', 'B*15', 'B*35', 'C*03', 'C*04'],
    'UDN08' :	['A*02', 'A*11', 'B*40', 'B*58', 'C*0202', 'C*07'],
    'UDN09' :	['A*02', 'A*68', 'B*14', 'B*27', 'C*0202', 'C*08'],
    'UDN10' :	['A*01', 'A*02', 'B*08', 'B*4402', 'C*07', 'C*0702'],
    'UDN11' :	['A*03', 'A*2402', 'B*35', 'B*45', 'C*04', 'C*16'],
    'UDN12' :	['A*2402', 'A*25', 'B*18', 'B*41', 'C*1202', 'C*17'],
    'UDN13' :	['A*01', 'A*25', 'B*13', 'B*39', 'C*0602', 'C*1202'],
    'UDN14' :	['A*03', 'A*29', 'B*4402', 'B*14', 'C*08', 'C*16'],
    'UDN15' :	['A*01', 'A*02', 'B*08', 'B*15', 'C*07', 'C*03'],
    'UDN16' :	['A*03', 'A*32', 'B*0702', 'B*27', 'C*07', 'C*0202'],
    'UDN17' :	['A*01', 'A*03', 'B*08', 'B*15', 'C*07', 'C*03'],
    'UDN18' :	['A*03', 'A*2402', 'B*0702', 'B*08', 'C*07'],
    'UDN19' :	['A*03', 'A*23', 'B*0702', 'C*07'],
    'UDN20' :	['A*02', 'A*31', 'B*27', 'B*4402', 'C*0202', 'C*05'],
    'UDN21' :	['A*03', 'B*08', 'B*35', 'C*07', 'C*04'],
    'UDN22' :	['A*03', 'B*0702', 'B*40', 'C*07', 'C*03'],
    'UDN23' :	['A*01', 'A*30', 'B*0702', 'B*18', 'C*07', 'C*05'],
    'UDN24' :	['A*01', 'A*32', 'B*08', 'B*13', 'C*07', 'C*0602'],
    'UDN25' :	['A*02', 'A*31', 'B*39', 'B*51', 'C*1202', 'C*1402'],
    'UDN26' :	['A*01', 'A*2402', 'B*52', 'B*57', 'C*1202', 'C*0602'],
    'UDN27' :	['A*01', 'A*03', 'B*08', 'B*14', 'C*07', 'C*08'],
    'UDN28' :	['A*26', 'A*68', 'B*4402', 'B*49', 'C*07', 'C*0702'],
    'UDN29' :	['A*02', 'A*11', 'B*18', 'B*58', 'C*07', 'C*0702'],
    'UDN30' :	['A*01', 'A*23', 'B*4402', 'B*49', 'C*04', 'C*07'],
    'UDN31' :	['A*02', 'A*2402', 'B*41', 'B*4402', 'C*17', 'C*05'],
    'UDN32' :	['A*01', 'A*2402', 'B*08', 'B*57', 'C*07'],
    'UDN33' :	['A*02', 'A*2402', 'B*18', 'B*49', 'C*07'],
    'UDN34' :	['A*01', 'A*02', 'B*08', 'B*4402', 'C*03', 'C*07'],
    'CLL_01' :	['A*01', 'A*02', 'B*08', 'B*27', 'C*07', 'C*0202'],
    'CLL_02' :	['A*02', 'B*15', 'B*56', 'C*03', 'C*0102'],
    'CLL_03' :	['A*11', 'B*18', 'B*39', 'C*05', 'C*1202'],
    'CLL_04' :  [],
    'CLL_05' :	['A*02', 'A*29', 'B*4402', 'C*05', 'C*16'],
    'CLL_06' :	['A*31', 'B*35', 'B*40', 'C*04', 'C*03'],
    'CLL_07' :	['A*02', 'A*2402', 'B*15', 'B*4402', 'C*03', 'C*16'],
    'CLL_08' :	['A*2402', 'A*25', 'B*18', 'B*49', 'C*07', 'C*1202'],
    'CLL_09' :	['A*03', 'A*2402', 'B*35', 'C*04'],
    'CLL_10' :	['A*02', 'A*2402', 'B*0702', 'B*4402', 'C*07', 'C*05'],
    'CLL_11' :	['A*01', 'A*26', 'B*08', 'B*38', 'C*07', 'C*1202'],
    'CLL_12' :	['A*01', 'A*02', 'B*13', 'B*4402', 'C*0602', 'C*16'],
    'CLL_13' :	['A*11', 'A*68', 'B*08', 'B*35', 'C*07', 'C*04'],
    'CLL_14' :	['A*01', 'A*2402', 'B*15', 'B*40', 'C*07', 'C*0202'],
    'CLL_15' :	['A*01', 'A*02', 'B*27', 'B*41', 'C*0202', 'C*17'],
    'CLL_16' :	['A*02', 'B*27', 'C*0202', 'C*1202'],
    'CLL_17' :	['A*01', 'A*11', 'B*08', 'B*52', 'C*07', 'C*1202'],
    'CLL_18' :	['A*03', 'B*0702', 'C*07'],
    'CLL_19' :	['A*03', 'A*26', 'B*0702', 'B*38', 'C*07', 'C*1202'],
    'CLL_20' :	['A*02', 'B*15', 'C*03', 'C*0602'],
    'CLL_21' :	['A*31', 'A*68', 'B*40', 'C*03'],
    'CLL_22' :	['A*02', 'A*2402', 'B*0702', 'B*13', 'C*07', 'C*0602'],
    'HNSCC' :	['A*2402', 'A*68', 'B*35', 'B*4402', 'C*04', 'C*07'],
    'JY':   	['A*0201', 'B*0702', 'C*0702'],
    'RCC':  	['A*11', 'A*2402', 'B*0702', 'B*51']
}

allele_map_II = { }

allele_map = allele_map_I

pin_file_dir = '/home/rw762/palmer_scratch/PXD038782/timstof/HLA_I_pin_mzml'
mzml_folder = Path('/home/rw762/palmer_scratch/PXD038782/timstof/HLA_I_pin_mzml')
output_folder = Path('/home/rw762/palmer_scratch/PXD038782/timstof/mhcbooster_I')
output_folder.mkdir(parents=True, exist_ok=True)

pin_files = list(Path(pin_file_dir).rglob('*.pin'))

rt_all_models = ['Deeplc_hela_hf', 'AlphaPeptDeep_rt_generic', 'Prosit_2019_irt', 'Prosit_2024_irt_cit']
ms2_all_models = ['Prosit_2019_intensity', 'Prosit_2024_intensity_cit', 'Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_CID', 'Prosit_2020_intensity_HCD',
                  'ms2pip_HCD2021', 'ms2pip_timsTOF2023', 'ms2pip_Immuno_HCD', 'ms2pip_timsTOF2024']
ccs_top_models = ['AlphaPeptDeep_ccs_generic']

for i in range(603,len(pin_files)):
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

    validator = MhcValidator(max_threads=40)  # Open a MHCvalidator instance, a new one has to be opened for each .pin file
    validator.set_mhc_params(alleles=alleles, mhc_class='I', min_pep_len=8)  # Load the alleles you specified above
    validator.load_data(pin, filetype='pin')  # Load the pin file
    validator.run(sequence_encoding=True, mhcflurry=use_mhc_scores, bigmhc=use_mhc_scores,
                  netmhcpan=use_mhc_scores, mixmhc2pred=use_mhc_scores,
                  autort=False, deeplc=False,
                  im2deep=True,
                  peptdeep=False,
                  # koina_predictors=rt_all_models + ms2_all_models + ccs_top_models,
                  koina_predictors=['Prosit_2024_irt_cit', 'Deeplc_hela_hf', 'AlphaPeptDeep_rt_generic',
                                    'ms2pip_timsTOF2024', 'Prosit_2023_intensity_timsTOF', 'AlphaPeptDeep_ccs_generic'],
                  fine_tune=False,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  report_directory=output_folder / f'{file_name}')
