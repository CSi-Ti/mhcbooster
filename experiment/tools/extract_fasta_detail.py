
import pandas as pd
from pathlib import Path
from src.utils.fasta_parser import parse_keys





if __name__ == '__main__':
    fasta_path = '/mnt/d/workspace/mhc-booster/data/JY_1M_0131/2025-01-29-decoys-contam-new_modified_fasta_GP2.fasta.fas'
    result_path = '/mnt/d/workspace/mhc-booster/experiment/JY_1M_0131/'
    key_info_map = parse_keys(fasta_path)
    for mhcb_path in Path(result_path).rglob('*.MhcValidator_annotated.tsv'):
        result_df = pd.read_csv(mhcb_path, sep='\t')
        prot_data = result_df['Proteins']
        result_df['UniqueIdentifier'] = None
        result_df['EntryName'] = None
        result_df['Description'] = None
        for i in range(len(prot_data)):
            protein_list = prot_data[i].split(';')
            description = set()
            uniqid = set()
            entryname = set()
            for protein in protein_list:
                protein = protein.replace('@', '').strip()
                if len(protein) == 0 or protein.startswith('rev_'):
                    continue
                prot_info = key_info_map[protein]
                description.add(prot_info)
                prot_split = prot_info.split('|')
                if len(prot_split) == 3:
                    info = prot_split[1]
                else:
                    info = prot_split[0]
                info_split = info.split('-')
                if len(info_split) == 3:
                    uniqid.add(info_split[1])
                    entryname.add(info_split[2])
                else:
                    if len(prot_split) == 3:
                        uniqid.add(prot_split[1])
                        entryname.add(prot_split[2])
            descriptions = ';'.join(description)
            uniqids = ';'.join(uniqid)
            entrynames = ';'.join(entryname)
            result_df.loc[i, 'UniqueIdentifier'] = uniqids
            result_df.loc[i, 'EntryName'] = entrynames
            result_df.loc[i, 'Description'] = descriptions
        output_path = mhcb_path.as_posix().replace('.tsv', '.ProtInfo.tsv')
        result_df.to_csv(output_path, sep='\t')