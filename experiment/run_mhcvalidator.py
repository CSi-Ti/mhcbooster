
import os
from pathlib import Path
from mhcvalidator.validator import MhcValidator

# alleles = ['HLA-A0201', 'HLA-B0702', 'HLA-C0702']
# pin_files = Path('/mnt/d/data/PXD052187/fragpipe/Search_EBV_0903').rglob('*.pin')
# output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/PXD052187/mhcvalidator_2')
# mzml_folder = Path('/mnt/d/data/PXD052187/mzml')

alleles = ['HLA-A0201', 'HLA-B0702', 'HLA-C0702']
pin_files = list(Path('/mnt/d/data/JY_1_10_25M/fragpipe/Search_NoBooster_FDR1').rglob('*.pin'))
output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/JY_1_10_25M/mhcflurry')
# mzml_folder = Path('/mnt/d/data/JY_1_10_25M/timsconvert')
mzml_folder = Path('/mnt/e/test')

# alleles = ['HLA-A0101', 'HLA-A2415', 'HLA-B5701', 'HLA-C0602']
# pin_files = Path('/mnt/d/data/HL-60/fragpipe/Search_human_msbooster').rglob('*.pin')
# output_folder = Path('/mnt/d/workspace/mhc-validator-2/experiment/HL-60/mhcvalidator_2_beta_all')
# mzml_folder = Path('/mnt/d/data/HL-60/raw')

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
# peptdeep_models = ['AlphaPeptDeep_ms2_generic'] # TimsTOF
# prosit_models = ['Prosit_2019_intensity', 'Prosit_2024_intensity_cit', 'Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_CID', 'Prosit_2020_intensity_HCD']
# ms2_bug_models = ['UniSpec', 'ms2pip_TTOF5600', 'Prosit_2024_intensity_PTMs_gl']
# ms2pip_models = ['ms2pip_HCD2021', 'ms2pip_timsTOF2023', 'ms2pip_iTRAQphospho', 'ms2pip_Immuno_HCD', 'ms2pip_timsTOF2024', 'ms2pip_CID_TMT']
# ms2pip_top2_models = ['ms2pip_timsTOF2024', 'ms2pip_HCD2021']
# prosit_top2_models = ['Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_HCD']
rt_top2_models = ['Prosit_2019_irt','Prosit_2024_irt_cit']
rt_all_models = ['Deeplc_hela_hf', 'AlphaPeptDeep_rt_generic', 'Prosit_2019_irt', 'Prosit_2024_irt_cit', 'Chronologer_RT']
ms2_all_models = ['Prosit_2019_intensity', 'Prosit_2024_intensity_cit', 'Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_CID', 'Prosit_2020_intensity_HCD',
                  'ms2pip_HCD2021', 'ms2pip_timsTOF2023', 'ms2pip_Immuno_HCD', 'ms2pip_timsTOF2024']
ccs_top_models = ['AlphaPeptDeep_ccs_generic']
ms2_top4_models = ['ms2pip_timsTOF2024', 'ms2pip_HCD2021', 'Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_HCD']
ms2_top2_tims_models = ['ms2pip_timsTOF2024', 'Prosit_2023_intensity_timsTOF']
ms2_top2_hcd_models = ['ms2pip_HCD2021', 'Prosit_2020_intensity_HCD']

for pin in pin_files:
    file_name = pin.stem
    if 'edited' in file_name:
        continue
    print(file_name)
    validator = MhcValidator(max_threads=os.cpu_count()//2) # Open a MHCvalidator instance, a new one has to be opened for each .pin file
    validator.set_mhc_params(alleles=alleles, mhc_class='I') # Load the alleles you specified above
    validator.load_data(pin, filetype='pin') # Load the pin file
    validator.run(sequence_encoding=False, netmhcpan=False, mhcflurry=True, bigmhc=False,
                  autort=False, deeplc=False,
                  im2deep=False,
                  peptdeep=False,
                  # koina_predictors=['Prosit_2019_irt','Prosit_2024_irt_cit', 'ms2pip_Immuno_HCD', 'Prosit_2020_intensity_HCD', 'Prosit_2024_intensity_cit'],
                  # koina_predictors=['IM2Deep','AlphaPeptDeep_ccs_generic'],
                  # koina_predictors=['UniSpec', 'ms2pip_HCD2021', 'ms2pip_timsTOF2023', 'ms2pip_iTRAQphospho', 'ms2pip_Immuno_HCD', 'ms2pip_TTOF5600', 'ms2pip_timsTOF2024', 'ms2pip_CID_TMT', 'AlphaPeptDeep_ms2_generic', 'Prosit_2019_intensity', 'Prosit_2024_intensity_cit', 'Prosit_2023_intensity_timsTOF', 'Prosit_2020_intensity_CID', 'Prosit_2024_intensity_XL_NMS2', 'Prosit_2023_intensity_XL_CMS2', 'Prosit_2020_intensity_TMT', 'Prosit_2020_intensity_HCD', 'Prosit_2023_intensity_XL_CMS3'],
                  # koina_predictors= rt_all_models + ccs_top_models + ms2_all_models,
                  koina_predictors= [],
                  # koina_predictors=['Prosit_2024_irt_cit', 'ms2pip_HCD2021', 'Prosit_2023_intensity_timsTOF'],
                  # koina_predictors=['Prosit_2024_irt_cit', 'ms2pip_HCD2021', 'Prosit_2019_intensity', 'Prosit_2020_intensity_HCD'],
                  # koina_predictors=['Prosit_2023_intensity_timsTOF'],
                  fine_tune=True,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  report_directory=output_folder / f'{file_name}_MhcValidator')
