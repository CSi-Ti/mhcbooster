

import pandas as pd
import numpy as np
from mhcbooster.utils.peptide import remove_previous_and_next_aa, remove_modifications, remove_charge
# pd.read_csv('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/human_yeast/mhcbooster/JY_Class1_25M_DDA_60min_Slot1-12_1_552_MHCBooster/peptide.tsv')

pin_path = '/mnt/d/data/JY_1_10_25M/fragpipe/Search_human_yeast/JY_Class1_25M_DDA_60min_Slot1_12_1_552/JY_Class1_25M_DDA_60min_Slot1-12_1_552.pin'

result_df = pd.read_csv(pin_path, sep='\t')
print(len(result_df))
result_df = result_df[-(result_df['Proteins'].str.contains('YEAST') * result_df['Proteins'].str.contains('HUMAN'))]
print(len(result_df))
result_df['Sequence'] = remove_modifications(remove_charge(remove_previous_and_next_aa(result_df['Peptide'].to_numpy())))
result_df = result_df[result_df['Label'] == 1]
print(np.sum(result_df['Proteins'].str.contains('HUMAN')))
print(np.sum(result_df['Proteins'].str.contains('YEAST')))
print(len(result_df[result_df['Proteins'].str.contains('HUMAN')]['Sequence'].unique()))
print(len(result_df[result_df['Proteins'].str.contains('YEAST')]['Sequence'].unique()))