
import os
from pathlib import Path
from src.main_mhcbooster import MhcValidator

alleles = ['HLA-A0201', 'HLA-B0702', 'HLA-C0702']
# alleles = ['DQB1*04', 'DQB1*06', 'DRB1*08', 'DRB1*13']
pin_files = Path('/mnt/d/data/JY_1_10_25M/fragpipe/Search_0226').rglob('*.pin')
mzml_folder = Path('/mnt/d/data/JY_1_10_25M/raw')
output_folder = Path('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/Search_0226')

# alleles = ['HLA-A0220', 'HLA-A6801', 'HLA-B3503', 'HLA-B3901', 'HLA-C0401', 'HLA-C0702']
# # alleles = ['HLA-A0201', 'HLA-B0702', 'HLA-C0702']
# pin_files = list(Path('/mnt/d/data/RA_Fractionation_Replicate_1/fragpipe/Search_0208/').rglob('*.pin'))
# pin_files = sorted(pin_files, key=lambda f: -f.stat().st_size)
# mzml_folder = Path('/mnt/d/data/RA_Fractionation_Replicate_1/raw')
# output_folder = Path('/mnt/d/workspace/mhc-booster/experiment/RA_Fractionation_Replicate_1/test')

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

auto_predict_predictor = False
rt_predictors = ['Prosit_2019_irt', 'Prosit_2024_irt_cit']
ms2_predictors = ['ms2pip_timsTOF2024', 'Prosit_2023_intensity_timsTOF']
ccs_predictors = ['AlphaPeptDeep_ccs_generic', 'IM2Deep']
app_predictors = ['mhcflurry', 'netmhcpan', 'bigmhc']
# app_predictors = ['netmhciipan', 'mixmhc2pred']
# app_predictors = []
for pin in pin_files:
    file_name = pin.stem
    if 'edited' in file_name:
        continue
    print(file_name)
    validator = MhcValidator(max_threads=os.cpu_count()//2) # Open a MHCvalidator instance, a new one has to be opened for each .pin file
    validator.set_mhc_params(alleles=alleles, mhc_class='I') # Load the alleles you specified above
    validator.load_data(pin, filetype='pin') # Load the pin file
    validator.run(sequence_encoding=True,
                  app_predictors=app_predictors,
                  auto_predict_predictor=auto_predict_predictor,
                  rt_predictors=rt_predictors,
                  ms2_predictors=ms2_predictors,
                  ccs_predictors=ccs_predictors,
                  fine_tune=False,
                  n_splits=5,
                  mzml_folder=mzml_folder,
                  fasta_path='/mnt/d/data/Library/2025-02-26-decoys-contam-JY_var_splicing_0226.fasta.fas',
                  report_directory=output_folder / f'{file_name}')

    if auto_predict_predictor:
        rt_predictors = validator.rt_predictors
        ms2_predictors = validator.ms2_predictors
        ccs_predictors = validator.ccs_predictors
        auto_predict_predictor = False

