
import os
import sys
from pathlib import Path

from mhcbooster.main_mhcbooster import MHCBooster, run_mhcbooster
from mhcbooster.report.combined_reporter import CombinedReporter

alleles = ['HLA-A0201', 'HLA-B0702', 'HLA-C0702'] # HLA-A0201 HLA-B0702 HLA-C0702
pin_files = Path('/mnt/d/data/JY_1_10_25M/fragpipe/Search_0226').rglob('*.pin')
mzml_folder = Path('/mnt/d/data/JY_1_10_25M/raw')
output_folder = Path('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/test')

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

auto_predict_predictor = False
fasta_path = '/mnt/d/data/JY_1_10_25M/2024-09-03-decoys-contam-Human.fasta'
rt_predictors = ['Prosit_2019_irt', 'Prosit_2024_irt_cit']
ms2_predictors = ['ms2pip_timsTOF2024', 'Prosit_2023_intensity_timsTOF']
ccs_predictors = ['IM2Deep']
app_predictors = ['mhcflurry', 'netmhcpan']



run_mhcbooster(pin_files, sequence_encoding=True, alleles=alleles, mhc_class='I', app_predictors=app_predictors,
    auto_predict_predictor=auto_predict_predictor, rt_predictors=rt_predictors, ms2_predictors=ms2_predictors,
    ccs_predictors=ccs_predictors, fine_tune=False, fasta_path=fasta_path, mzml_folder=mzml_folder,
    output_folder=output_folder)