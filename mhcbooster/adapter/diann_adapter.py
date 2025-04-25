
import re
import numpy as np
import pandas as pd
from pathlib import Path

from mhcbooster.utils.constants import UNIMOD_MASS_MAP, PROTON_MASS


def convert_to_pin(file_path, output_path):
    file_path = Path(file_path)
    if file_path.suffix == '.tsv':
        result_df = pd.read_csv(file_path, sep='\t')
    elif file_path.suffix == '.parquet':
        result_df = pd.read_parquet(file_path)
    else:
        raise ValueError(f'Unsupported DIA-NN result file type: {file_path.name}')

    # replace unimod to mass diff
    def convert_unimod_to_mass(seq):
        matches = list(re.finditer(r"\(UniMod:(\d+)\)", seq))
        for m in reversed(matches):
            unimod_id = m.group(1)
            if unimod_id in UNIMOD_MASS_MAP.keys():
                if m.start() == 0:
                    seq = 'n[' + str(UNIMOD_MASS_MAP[unimod_id]) + ']' + seq[m.end():]
                else:
                    seq = seq[:m.start()] + '[' + str(UNIMOD_MASS_MAP[unimod_id]) + ']' + seq[m.end():]
            else:
                print(f'Ignored supported Unimod: {unimod_id}.')
                seq = seq[:m.start()] + seq[m.end():]
        return seq
    result_df['SpecId'] = result_df['Precursor.Id']
    result_df['Peptide'] = result_df['Modified.Sequence'].apply(convert_unimod_to_mass)
    result_df['Peptide'] = '-.' + result_df['Peptide'] + result_df['Precursor.Charge'].astype(str) + '.-'
    result_df['Label'] = 1 - 2 * result_df['Decoy'] #TODO percolator
    result_df['Proteins'] = result_df['Protein.Ids']
    result_df.loc[result_df['Decoy'] == 1, 'Proteins'] = 'rev_' + result_df.loc[result_df['Decoy'] == 1, 'Proteins']
    result_df['Charge'] = result_df['Precursor.Charge']
    result_df['ExpMass'] = (result_df['Precursor.Mz'] - PROTON_MASS) * result_df['Charge']
    result_df['qvalue'] = result_df['Q.Value']
    result_df['ScanNr'] = np.arange(1, len(result_df) + 1)
    result_df = result_df.drop(columns=['Modified.Sequence', 'Decoy', 'Protein.Ids', 'Protein.Group', 'Protein.Sites',
                                        'Protein.Names', 'Genes', 'Precursor.Charge', 'Site.Occupancy.Probabilities',
                                        'Run.Index', 'Run', 'Channel', 'Precursor.Id', 'Stripped.Sequence'])
    header_cols = ['SpecId', 'Label', 'ScanNr', 'ExpMass']
    result_df = result_df[header_cols + [c for c in result_df.columns if c not in header_cols]]
    result_df.to_csv(output_path, sep='\t', index=False)


if __name__ == '__main__':
    file_path = '/mnt/d/data/JY_DIA/test/report.parquet'
    output_path = '/mnt/d/data/JY_DIA/test/report.pin'
    convert_to_pin(file_path, output_path)