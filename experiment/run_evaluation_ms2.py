import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyteomics.parser import length
from src.utils.peptide import replace_uncommon_aas, remove_charge, remove_previous_and_next_aa, remove_modifications
from src.predictors.netmhcpan_helper import NetMHCpanHelper

N_THREADS = os.cpu_count() // 2


def read_feature(feature_path):

    feature_df = pd.read_csv(feature_path, sep='\t')

    pep_df = feature_df[feature_df['mhcv_label'] == 1]
    pep_df = pep_df[pep_df['mhcv_q-value'] <= 0.01]
    # entropy_score = pep_df['ms2pip_timsTOF2024_entropy_score']
    # for score in entropy_score:
    #     print(score)
    el_columns = [column for column in feature_df.columns if column.endswith('ELScore') and 'log' not in column]
    for i in range(len(pep_df)):
        max_el = max(pep_df.iloc[i][el_columns])
        print(max_el)

if __name__ == '__main__':

    read_feature('/mnt/d/workspace/mhc-booster/experiment/MSV000091456/A375_lowInput_IP/HLA-I/mhcbooster/E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep4_Slot2-5_1_3522/E_20221201_NO30_400nL_HLAc1_4e7_directIP_titration_rep4_Slot2-5_1_3522.features.tsv')
