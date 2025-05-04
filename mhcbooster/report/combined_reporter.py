

import tempfile
import subprocess
import numpy as np
import pandas as pd

from pathlib import Path
from copy import deepcopy
from pyteomics import fasta, protxml

class CombinedReporter:

    def __init__(self, result_folder, fasta_path, pep_fdr=0.01, seq_fdr=0.01, infer_protein=True,
                 decoy_prefix = 'rev_', remove_contaminant=False):

        self.result_folder = Path(result_folder)
        self.fasta_path = fasta_path
        self.pep_fdr = pep_fdr
        self.seq_fdr = seq_fdr
        self.infer_protein = infer_protein
        self.decoy_prefix = decoy_prefix
        self.remove_contaminant = remove_contaminant

        if self.fasta_path is None:
            fasta_paths = list(self.result_folder.glob('*.fasta')) + list(self.result_folder.glob('*.fas'))
            if len(fasta_paths) > 0:
                self.fasta_path = fasta_paths[0]

        self.fasta_map = None
        self.seq_prot_df = None
        self.combined_high_confidence_sequences = None

    def prepare_fasta_map(self):
        # Prepare fasta: some sequences are missing in prot.xml
        fasta_map = {}
        if self.fasta_path is not None:
            for protein in fasta.read(self.fasta_path):
                description = protein.description.strip()
                protein_name = description.split(' ')[0]
                protein_description = description[len(protein_name) + 1:]
                fasta_map[protein_name] = protein_description
        return fasta_map


    def do_protein_inference(self):
        if self.fasta_path is None:
            return pd.DataFrame()

        ### Run ProteinProphet
        pep_xml_paths = [path.absolute().as_posix() for path in self.result_folder.rglob('peptide.pep.xml')]
        with tempfile.NamedTemporaryFile('w') as pep_xml_list:
            pep_xml_list.write('\n'.join(pep_xml_paths))
            pep_xml_list.flush()
            philosopher_exe_path = Path(__file__).parent.parent / 'third_party' / 'philosopher_v5.1.0_linux_amd64' / 'philosopher'
            subprocess.run(f'{philosopher_exe_path} workspace --init', cwd=self.result_folder, shell=True)
            subprocess.run(f'{philosopher_exe_path} proteinprophet --maxppmdiff 2000000 --output combined {pep_xml_list.name}', cwd=self.result_folder, shell=True)
            subprocess.run(f'{philosopher_exe_path} workspace --clean --nocheck', cwd=self.result_folder, shell=True)


        # Generate sequence-protein map from prot.xml
        if not (self.result_folder / 'combined.prot.xml').exists():
            print('Protein inference is skipped.')
            return pd.DataFrame()
        prot_data = list(protxml.read(str(self.result_folder / 'combined.prot.xml')))
        seq_prot_map = {}
        for prot in prot_data:
            protein = prot['protein'][0]
            prot_desc_split = [t for t in protein['protein_description'].split(' ') if t.startswith('GN=')]
            gene_name = 'UNANNOTATED' if len(prot_desc_split) == 0 else prot_desc_split[0].replace('GN=', '')

            protein_list = [{
                'protein_name': protein['protein_name'],
                'protein_description': protein['protein_description'],
                'gene_name': self.decoy_prefix + gene_name if protein['protein_name'].startswith(
                    self.decoy_prefix) else gene_name,
                'n_related_peptides': len(protein['peptide'])
            }]
            if protein['n_indistinguishable_proteins'] > 1:
                for dup_prot in protein['indistinguishable_protein']:
                    prot_desc_split = [t for t in dup_prot['protein_description'].split(' ') if t.startswith('GN=')]
                    gene_name = 'UNANNOTATED' if len(prot_desc_split) == 0 else prot_desc_split[0].replace('GN=', '')
                    protein_list.append({
                        'protein_name': dup_prot['protein_name'],
                        'protein_description': dup_prot['protein_description'],
                        'gene_name': self.decoy_prefix + gene_name if dup_prot['protein_name'].startswith(
                            self.decoy_prefix) else gene_name,
                        'n_related_peptides': len(protein['peptide'])
                    })
            for peptide in prot['protein'][0]['peptide']:
                sequence = peptide['peptide_sequence']
                group_weight = peptide['group_weight']
                n_sibling_peptides = peptide['n_sibling_peptides']
                for prot_map in protein_list:
                    prot_map['group_weight'] = group_weight
                    prot_map['n_sibling_peptides'] = n_sibling_peptides

                if sequence not in seq_prot_map.keys():
                    seq_prot_map[sequence] = []
                protein_names = [protein['protein_name'] for protein in seq_prot_map[sequence]]
                for protein in protein_list:
                    if protein['protein_name'] not in protein_names:
                        seq_prot_map[sequence].append(deepcopy(protein))
                        protein_names.append(protein['protein_name'])

        # find the best protein for each sequence
        seq_idx_map = {}
        for sequence in seq_prot_map.keys():
            protein_list = seq_prot_map[sequence]
            max_weight = max(protein['group_weight'] for protein in protein_list)
            indices = [i for i, protein in enumerate(protein_list) if protein['group_weight'] == max_weight]
            if len(indices) == 1:
                seq_idx_map[sequence] = indices[0]
                continue

            n_related_peptides = [protein_list[i]['n_related_peptides'] for i in indices]
            max_related_peptides = max(n_related_peptides)
            indices = [indices[i] for i in range(len(n_related_peptides)) if
                       n_related_peptides[i] == max_related_peptides]
            if len(indices) == 1:
                seq_idx_map[sequence] = indices[0]
                continue

            protein_names = [protein_list[i]['protein_name'] for i in indices]
            index = indices[protein_names.index(min(protein_names))]
            seq_idx_map[sequence] = index

        sequences = list(seq_prot_map.keys())
        proteins = [None] * len(sequences)
        protein_ids = [None] * len(sequences)
        entry_names = [None] * len(sequences)
        genes = [None] * len(sequences)
        protein_descriptions = [None] * len(sequences)
        mapped_proteins = [None] * len(sequences)
        mapped_genes = [None] * len(sequences)

        # Use list comprehensions and set operations to avoid redundant calculations
        for i, sequence in enumerate(sequences):
            protein_list = seq_prot_map[sequence]
            best_index = seq_idx_map[sequence]
            best_protein = protein_list[best_index]

            protein_name = best_protein['protein_name']
            proteins[i] = protein_name
            protein_split = protein_name.split('|')
            protein_ids[i] = self.decoy_prefix + protein_split[1] if protein_name.startswith(self.decoy_prefix) else protein_split[1]
            entry_names[i] = self.decoy_prefix + protein_split[2] if protein_name.startswith(self.decoy_prefix) else protein_split[2]
            genes[i] = best_protein['gene_name']
            protein_descriptions[i] = best_protein['protein_description']

            mapped_proteins_list = list(set([protein_list[j]['protein_name'] for j in range(len(protein_list)) if j != best_index]))
            mapped_target_proteins = [protein_name for protein_name in mapped_proteins_list if not protein_name.startswith(self.decoy_prefix)]
            mapped_decoy_proteins = [protein_name for protein_name in mapped_proteins_list if protein_name.startswith(self.decoy_prefix)]
            mapped_proteins[i] = ','.join(sorted(mapped_target_proteins) + sorted(mapped_decoy_proteins))

            mapped_genes_list = list(set([protein_list[j]['gene_name'] for j in range(len(protein_list)) if j != best_index]))
            mapped_target_genes = [gene_name for gene_name in mapped_genes_list if not gene_name.startswith(self.decoy_prefix)]
            mapped_decoy_genes = [gene_name for gene_name in mapped_genes_list if gene_name.startswith(self.decoy_prefix)]
            mapped_genes[i] = ','.join(sorted(mapped_target_genes) + sorted(mapped_decoy_genes))

        seq_prot_df = pd.DataFrame({
            'sequence': sequences,
            'protein': proteins,
            'protein_id': protein_ids,
            'entry_name': entry_names,
            'gene': genes,
            'protein_description': protein_descriptions,
            'mapped_protein': mapped_proteins,
            'mapped_gene': mapped_genes
        })
        return seq_prot_df


    def update_proteins(self, file_name):
        if self.fasta_path is None:
            return
        for result_path in self.result_folder.rglob(file_name):
            print(f'Inserting protein inference results to {result_path.name}')
            result_df = pd.read_csv(result_path, sep='\t')
            protein_idx = result_df.columns.get_loc('protein')
            protein_col = result_df['protein']
            result_df = result_df.drop(['protein', 'protein_id', 'entry_name', 'gene', 'protein_description',
                                        'mapped_protein', 'mapped_gene'], axis=1)
            merged_df = pd.DataFrame()
            merged_df['sequence'] = result_df['sequence']
            merged_df = pd.merge(merged_df, self.seq_prot_df, how='left', on='sequence')
            for i in range(1, len(merged_df.columns)):
                result_df.insert(protein_idx + i - 1, merged_df.columns[i], merged_df.iloc[:, i])

            for i, row in result_df.iterrows():
                if pd.isna(row['protein']):
                    proteins = sorted([protein for protein in protein_col[i].split(';') if len(protein.strip()) > 0])
                    if len(proteins) == 0:
                        continue
                    protein = proteins[0]
                    protein_description = ''
                    gene = ''
                    protein_id, entry_name = '', ''
                    if protein in self.fasta_map.keys():
                        protein_description = self.fasta_map[protein]
                        prot_desc_split = [t for t in protein_description.split(' ') if t.startswith('GN=')]
                        gene = 'UNANNOTATED' if len(prot_desc_split) == 0 else prot_desc_split[0][3:]
                        gene = self.decoy_prefix + gene if protein.startswith(self.decoy_prefix) else gene
                        prot_name_split = protein.split('|')
                        if len(prot_name_split) == 3:
                            protein_id = self.decoy_prefix + prot_name_split[1] if protein.startswith(self.decoy_prefix) else prot_name_split[1]
                            entry_name = self.decoy_prefix + prot_name_split[2] if protein.startswith(self.decoy_prefix) else prot_name_split[2]

                    mapped_proteins = proteins[1:]
                    mapped_target_proteins = [protein_name for protein_name in mapped_proteins if not protein_name.startswith(self.decoy_prefix)]
                    mapped_decoy_proteins = [protein_name for protein_name in mapped_proteins if protein_name.startswith(self.decoy_prefix)]
                    mapped_proteins = sorted(mapped_target_proteins) + sorted(mapped_decoy_proteins)
                    mapped_genes = []
                    for mapped_protein in mapped_proteins:
                        description = self.fasta_map[mapped_protein]
                        prot_desc_split = [t for t in description.split(' ') if t.startswith('GN=')]
                        gene_name = 'UNANNOTATED' if len(prot_desc_split) == 0 else prot_desc_split[0][3:]
                        gene_name = self.decoy_prefix + gene_name if mapped_protein.startswith(self.decoy_prefix) else gene_name
                        mapped_genes.append(gene_name)
                    mapped_genes = list(set(mapped_genes))
                    mapped_protein = ','.join(mapped_proteins)
                    mapped_gene = ','.join(mapped_genes)
                    result_df.loc[i, ['protein', 'protein_description', 'gene', 'protein_id', 'entry_name', 'mapped_protein', 'mapped_gene']] = \
                        [protein, protein_description, gene, protein_id, entry_name, mapped_protein, mapped_gene]
            result_df.to_csv(result_path, sep='\t', index=False)


    def remove_contaminants(self, file_name):
        contam_fasta = fasta.read(str(Path(__file__).parent / 'contaminant.fasta.fas'))
        contam_proteins = [p.description.split(' ')[0] for p in contam_fasta]
        for result_path in self.result_folder.rglob(file_name):
            print(f'Removing contaminants in {result_path.parent.name} {result_path.name}', end='\t')
            result_df = pd.read_csv(result_path, sep='\t')
            result_df_wo_contam = result_df[-result_df['protein'].isin(contam_proteins)]
            print(f'{len(result_df) - len(result_df_wo_contam)} {file_name.split(".")[0]}s removed.')
            result_df_wo_contam.to_csv(result_path, sep='\t', index=False)


    def combine_result(self, file_name):

        common_cols = ['peptide', 'sequence', 'prev_AA', 'next_AA', 'seq_len', 'charge',
                       'protein', 'protein_id', 'entry_name', 'gene', 'protein_description',
                       'mapped_protein', 'mapped_gene']
        if 'peptide' in file_name:
            fdr = self.pep_fdr
            group_key = 'peptide'
        else:
            fdr = self.seq_fdr
            group_key = 'sequence'
        print(f'Generating combined {group_key} result to combined_{group_key}.tsv')

        result_paths = list(self.result_folder.rglob(file_name))
        all_dfs = []
        for result_path in result_paths:
            file_name = Path(result_path).parent.name
            result_df = pd.read_csv(result_path, sep='\t')
            fdr_col = [col for col in result_df.columns if '_qvalue' in col][0]
            result_df = result_df[(result_df['label'] == 'Target') * (result_df[fdr_col] <= fdr)]
            if not pd.api.types.is_string_dtype(result_df['charge']):
                result_df['charge'] = result_df['charge'].astype(str)

            cols_to_drop = set(result_df.columns) & {'label', 'score', 'min_rank', fdr_col}
            cols_to_drop |= {col for col in result_df.columns if '_binder' in col}
            result_df.drop(columns=cols_to_drop, inplace=True, errors='ignore')
            result_df = result_df.rename(columns={'spectral_count': file_name + '_spectral_count'})
            if 'binder' in result_df.columns:
                result_df = result_df.rename(columns={'binder': file_name + '_binder'})
                result_df = result_df.rename(columns={'best_allele': file_name + '_allele'})
            all_dfs.append(result_df)

        combined_df = pd.concat(all_dfs, ignore_index=True).merge(pd.DataFrame(columns=common_cols), on=common_cols, how='outer')
        agg_dict = {
            'peptide': lambda x: ','.join(sorted(set(','.join(x.dropna()).split(',')))),
            'sequence': 'first',
            'prev_AA': lambda x: ','.join(sorted(set(','.join(x.dropna()).split(',')))),
            'next_AA': lambda x: ','.join(sorted(set(','.join(x.dropna()).split(',')))),
            'seq_len': 'first',
            'charge': lambda x: ','.join(sorted(set(','.join(x.dropna()).split(',')))),
            'protein': 'first',
            'protein_id': 'first',
            'entry_name': 'first',
            'gene': 'first',
            'protein_description': 'first',
            'mapped_protein': 'first',
            'mapped_gene': 'first'
        }
        agg_dict.pop(group_key)
        for col in combined_df.columns[13:]:
            agg_dict[col] = lambda x: x.dropna().iloc[0] if not x.dropna().empty else np.nan
        combined_df = combined_df.groupby(group_key, as_index=False).agg(agg_dict)
        binder_cols = [col for col in combined_df.columns if 'binder' in col]
        allele_cols = [col for col in combined_df.columns if 'allele' in col]
        spectral_count_cols = [col for col in combined_df.columns if 'spectral_count' in col]
        for col in spectral_count_cols:
            combined_df[col] = combined_df[col].astype('Int64')
        cols = list(combined_df.columns[:13]) + spectral_count_cols + binder_cols + allele_cols
        combined_df = combined_df[cols]

        combined_len = len(combined_df)
        combined_df = combined_df[combined_df['sequence'].isin(self.combined_high_confidence_sequences)]
        print(f'Saving {len(combined_df)} {group_key}s to combined_{group_key}.tsv ({combined_len} before filtering)')
        combined_df.to_csv(self.result_folder / f'combined_{group_key}.tsv', sep='\t', index=False)
        print('Done.')


    def get_philosopher_reference(self):
        peptide_paths = list(self.result_folder.rglob('peptide.pep.xml'))
        if self.fasta_path is None or len(peptide_paths) == 0:
            return
        philosopher_exe_path = Path(__file__).parent.parent / 'third_party' / 'philosopher_v5.1.0_linux_amd64' / 'philosopher'
        for peptide_path in peptide_paths:
            sample_path = peptide_path.parent
            subprocess.run(f'{philosopher_exe_path} workspace --init', cwd=sample_path, shell=True)
            subprocess.run(f'{philosopher_exe_path} database --annotate {self.fasta_path}', cwd=sample_path, shell=True)
            subprocess.run(f'{philosopher_exe_path} filter --sequential --prot 1 --pep {self.pep_fdr} --tag rev_ --pepxml peptide.pep.xml --protxml ../combined.prot.xml --razor', cwd=sample_path, shell=True)

        pepxml_paths = ' '.join([str(path) for path in peptide_paths])
        sample_names = ' '.join([path.parent.name for path in peptide_paths])
        subprocess.run(f'{philosopher_exe_path} workspace --init', cwd=self.result_folder, shell=True)
        subprocess.run(f'{philosopher_exe_path} iprophet --decoy rev_ --nonsp --output combined {pepxml_paths}', cwd=self.result_folder, shell=True)
        subprocess.run(f'{philosopher_exe_path} abacus --razor --reprint --tag rev_ --protein --peptide {sample_names}', cwd=self.result_folder, shell=True)

        for peptide_path in peptide_paths:
            sample_path = peptide_path.parent
            subprocess.run(f'{philosopher_exe_path} workspace --clean --nocheck', cwd=sample_path, shell=True)
        subprocess.run(f'{philosopher_exe_path} workspace --clean --nocheck', cwd=self.result_folder, shell=True)
        try:
            reference_df = pd.read_csv(self.result_folder / 'combined_peptide.tsv', sep='\t')
            if 'Sequence' in reference_df.columns:
                self.combined_high_confidence_sequences = reference_df['Sequence'].to_numpy()
        except FileNotFoundError:
            print('Failed to generate the combined peptide list using Philosopher.')
            return


    def run(self):

        if self.infer_protein:
            self.fasta_map = self.prepare_fasta_map()
            self.seq_prot_df = self.do_protein_inference()
            self.update_proteins('psm.tsv')
            self.update_proteins('peptide.tsv')
            self.update_proteins('sequence.tsv')

        if self.remove_contaminant:
            self.remove_contaminants('psm.tsv')
            self.remove_contaminants('peptide.tsv')
            self.remove_contaminants('sequence.tsv')

        self.get_philosopher_reference()
        self.combine_result('peptide.tsv')
        self.combine_result('sequence.tsv')


if __name__ == '__main__':
    combined_reporter = CombinedReporter(result_folder='/mnt/d/data/JY_Fractionation_Replicate_1/mhcbooster_0307_copy',
                                         fasta_path='/mnt/d/data/JY_1_10_25M/2024-09-03-decoys-contam-Human_EBV_GD1_B95.fasta',
                                         infer_protein=False,
                                         remove_contaminant=True)
    combined_reporter.run()