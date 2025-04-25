from pathlib import Path

import pandas as pd

from mhcbooster.utils.constants import PROTON_MASS


def convert_result_to_pin(result_path):
    result_df = pd.read_csv(result_path, sep='\t')
    result_df = result_df.rename(columns={
        'Scan_No': 'ScanNr',
        'Sequence': 'Peptide',
        'File_Name': 'SpecId',
        'Exp.MH+': 'ExpMass'
    })
    result_df = result_df[-pd.isna(result_df['Target/Decoy'])]
    result_df = result_df[result_df['Q-value'] < 0.5]
    result_df.loc[result_df['Target/Decoy'] == 'target', 'Label'] = 1
    result_df.loc[result_df['Target/Decoy'] == 'decoy', 'Label'] = -1
    result_df['ExpMass'] -= PROTON_MASS
    result_df.drop(columns=['Modification', 'Positions', 'Target/Decoy', 'Others'], inplace=True)
    result_df.to_csv(result_path.parent/'JY_Class1_1M_DDA_60min_Slot1-10_1_541.pin', sep='\t', index=False) # TODO

if __name__ == '__main__':
    result_path = Path('/mnt/d/pFindWorkspace/single_1M/result/pFind.spectra')
    convert_result_to_pin(result_path)