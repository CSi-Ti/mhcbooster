import subprocess

import numpy as np
import pandas as pd
import matplotlib.backends.backend_pdf as plt_pdf

from pathlib import Path
from copy import deepcopy
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.cm import get_cmap
from pyteomics.mass import calculate_mass
from pyteomics import protxml, fasta
from src.utils.fdr import calculate_qs, calculate_peptide_level_qs, calculate_roc
from mhcnames import normalize_allele_name

class RunReporter:

    def __init__(self, report_directory, file_name, decoy_prefix):

        self.report_directory = Path(report_directory)
        self.file_name = file_name
        self.decoy_prefix = decoy_prefix

        self.psm_df = pd.DataFrame()
        self.peptide_df = pd.DataFrame()
        self.sequence_df = pd.DataFrame()


    def add_run_result(self, peptides, sequences, prev_aas, next_aas, labels, charges, scores, proteins):
        psm_qvalues = calculate_qs(scores, labels)
        pep_qvalues, _, _, peps, _ = calculate_peptide_level_qs(scores, labels, peptides)
        seq_qvalues, _, _, seqs, _ = calculate_peptide_level_qs(scores, labels, sequences)
        pep_qvalue_lookup = {pep: q for pep, q in zip(peps, pep_qvalues)}
        seq_qvalue_lookup = {seq: q for seq, q in zip(seqs, seq_qvalues)}

        self.psm_df['peptide'] = peptides
        self.psm_df['sequence'] = sequences
        self.psm_df['prev_AA'] = prev_aas
        self.psm_df['next_AA'] = next_aas
        self.psm_df['label'] = np.array(['Target' if label == 1 else 'Decoy' for label in labels])
        self.psm_df['seq_len'] = self.psm_df['sequence'].str.len().astype(str)
        self.psm_df['charge'] = charges
        self.psm_df['score'] = scores
        self.psm_df['psm_qvalue'] = psm_qvalues
        self.psm_df['pep_qvalue'] = np.array([pep_qvalue_lookup[pep] for pep in peptides])
        self.psm_df['seq_qvalue'] = np.array([seq_qvalue_lookup[seq] for seq in sequences])
        self.psm_df['protein'] = proteins
        self.psm_df['protein_id'] = ''
        self.psm_df['entry_name'] = ''
        self.psm_df['gene'] = ''
        self.psm_df['protein_description'] = ''
        self.psm_df['mapped_protein'] = ''
        self.psm_df['mapped_gene'] = ''


    def infer_protein(self, fasta_path, score_threshold=0):
        if fasta_path is None:
            return
        psm_df = self.psm_df[self.psm_df['score'] >= score_threshold]
        fasta_map = {}
        for protein in fasta.read(fasta_path):
            description = protein.description.strip()
            protein_name = description.split(' ')[0]
            protein_description = description[len(protein_name) + 1:]
            fasta_map[protein_name] = protein_description

        header = ['<?xml version="1.0" encoding="UTF-8"?>\n',
                  '<?xml-stylesheet type="text/xsl" href="pepXML_std.xsl"?>\n',
                  '<msms_pipeline_analysis xmlns="http://regis-web.systemsbiology.net/pepXML" xsi:schemaLocation="http://sashimi.sourceforge.net/schema_revision/pepXML/pepXML_v122.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n',
                  '<analysis_summary analysis="MHCBooster">\n',
                  '</analysis_summary>\n',
                  '<msms_run_summary>\n']

        with open(self.report_directory / 'peptide.pep.xml', 'w') as pepxml:
            pepxml.writelines(header)

            # fasta
            pepxml.write('<search_summary>\n')
            pepxml.write(f'<search_database local_path="{fasta_path}" type="AA"/>\n')
            pepxml.write(f'</search_summary>\n')

            for i, psm in psm_df.iterrows():
                sequence = psm['sequence']
                charge = psm['charge']
                proteins = [protein for protein in psm['protein'].split(';') if len(protein.strip()) > 0]
                score = psm['score']
                # score = 1 - psm['psm_qvalue']

                pepxml.write(f'<spectrum_query assumed_charge="{charge}" spectrum="{i}">\n')
                pepxml.write('<search_result>\n')
                pepxml.write(f'<search_hit peptide="{sequence}" calc_neutral_pep_mass="{calculate_mass(sequence)}" num_tot_proteins="{len(proteins)}" protein="{proteins[0]}">\n')
                for i in range(1, len(proteins)):
                    pepxml.write(f'<alternative_protein protein="{proteins[i]}"/>\n')
                pepxml.write('<analysis_result analysis="peptideprophet">\n')
                pepxml.write(f'<peptideprophet_result probability="{score}" all_ntt_prob="({score},{score},{score})">\n')
                pepxml.write('</peptideprophet_result>\n')
                pepxml.write('</analysis_result>\n')
                pepxml.write('</search_hit>\n')
                pepxml.write('</search_result>\n')
                pepxml.write('</spectrum_query>\n')

            pepxml.write('</msms_run_summary>\n')
            pepxml.write('</msms_pipeline_analysis>\n')

        philosopher_exe_path = Path(__file__).parent.parent.parent / 'third_party' / 'philosopher_v5.1.0_linux_amd64' / 'philosopher'
        subprocess.run(f'{philosopher_exe_path} workspace --init', cwd=self.report_directory, shell=True)
        subprocess.run(f'{philosopher_exe_path} proteinprophet --maxppmdiff 2000000 peptide.pep.xml', cwd=self.report_directory, shell=True)

        prot_data = list(protxml.read(str(self.report_directory / 'interact.prot.xml')))
        seq_prot_map = {}
        for prot in prot_data:
            protein = prot['protein'][0]
            prot_desc_split = [t for t in protein['protein_description'].split(' ') if t.startswith('GN=')]
            gene_name = 'UNANNOTATED' if len(prot_desc_split) == 0 else prot_desc_split[0].replace('GN=', '')

            protein_list = [{
                'protein_name': protein['protein_name'],
                'protein_description': protein['protein_description'],
                'gene_name': self.decoy_prefix + gene_name if protein['protein_name'].startswith(self.decoy_prefix) else gene_name,
                'n_related_peptides': len(protein['peptide'])
            }]
            if protein['n_indistinguishable_proteins'] > 1:
                for dup_prot in protein['indistinguishable_protein']:
                    prot_desc_split = [t for t in dup_prot['protein_description'].split(' ') if t.startswith('GN=')]
                    gene_name = 'UNANNOTATED' if len(prot_desc_split) == 0 else prot_desc_split[0].replace('GN=', '')
                    protein_list.append({
                        'protein_name': dup_prot['protein_name'],
                        'protein_description': dup_prot['protein_description'],
                        'gene_name': self.decoy_prefix + gene_name if dup_prot['protein_name'].startswith(self.decoy_prefix) else gene_name,
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
            indices = [indices[i] for i in range(len(n_related_peptides)) if n_related_peptides[i] == max_related_peptides]
            if len(indices) == 1:
                seq_idx_map[sequence] = indices[0]
                continue

            protein_names = [protein_list[i]['protein_name'] for i in indices]
            index = indices[protein_names.index(min(protein_names))]
            seq_idx_map[sequence] = index

        # prepare seq_prot_map for sequence-protein mapping
        for sequence in seq_prot_map.keys():
            protein_list = seq_prot_map[sequence]
            best_index = seq_idx_map[sequence]
            best_protein = protein_list[best_index]

            mapped_proteins = [protein_list[i]['protein_name'] for i in range(len(protein_list)) if i != best_index]
            mapped_proteins = list(set(mapped_proteins))
            mapped_target_proteins = [protein_name for protein_name in mapped_proteins if not protein_name.startswith(self.decoy_prefix)]
            mapped_decoy_proteins = [protein_name for protein_name in mapped_proteins if protein_name.startswith(self.decoy_prefix)]
            mapped_proteins = sorted(mapped_target_proteins) + sorted(mapped_decoy_proteins)

            mapped_genes = [protein_list[i]['gene_name'] for i in range(len(protein_list)) if i != best_index]
            mapped_genes = list(set(mapped_genes))
            mapped_target_genes = [gene_name for gene_name in mapped_genes if not gene_name.startswith(self.decoy_prefix)]
            mapped_decoy_genes = [gene_name for gene_name in mapped_genes if gene_name.startswith(self.decoy_prefix)]
            mapped_genes = sorted(mapped_target_genes) + sorted(mapped_decoy_genes)

            seq_prot_map[sequence] = {'best_protein': best_protein, 'mapped_proteins': mapped_proteins, 'mapped_genes': mapped_genes}

        for i, psm in self.psm_df.iterrows():
            sequence = psm['sequence']
            if sequence not in seq_prot_map.keys():
                proteins = sorted([protein for protein in psm['protein'].split(';') if len(protein.strip()) > 0])
                if len(proteins) == 0:
                    continue
                best_protein_name = proteins[0]
                best_description = fasta_map[best_protein_name]
                prot_desc_split = [t for t in best_description.split(' ') if t.startswith('GN=')]
                best_gene_name = 'UNANNOTATED' if len(prot_desc_split) == 0 else prot_desc_split[0].replace('GN=', '')
                mapped_proteins = proteins[1:]
                mapped_genes = []
                for mapped_protein in mapped_proteins:
                    description = fasta_map[mapped_protein]
                    prot_desc_split = [t for t in description.split(' ') if t.startswith('GN=')]
                    gene_name = 'UNANNOTATED' if len(prot_desc_split) == 0 else prot_desc_split[0].replace('GN=', '')
                    if mapped_protein.startswith(self.decoy_prefix):
                        gene_name = self.decoy_prefix + gene_name
                    mapped_genes.append(gene_name)
                mapped_genes = list(set(mapped_genes))
            else:
                best_protein_name = seq_prot_map[sequence]['best_protein']['protein_name']
                best_description = seq_prot_map[sequence]['best_protein']['protein_description']
                best_gene_name = seq_prot_map[sequence]['best_protein']['gene_name']
                mapped_proteins = seq_prot_map[sequence]['mapped_proteins']
                mapped_genes = seq_prot_map[sequence]['mapped_genes']

            self.psm_df.loc[i, 'protein'] = best_protein_name
            is_decoy = best_protein_name.startswith(self.decoy_prefix)
            protein_split = best_protein_name.split('|')
            if len(protein_split) == 3:
                protein_id = protein_split[1]
                entry_name = protein_split[2]
                self.psm_df.loc[i, 'protein_id'] = (self.decoy_prefix + protein_id) if is_decoy else protein_id
                self.psm_df.loc[i, 'entry_name'] = (self.decoy_prefix + entry_name) if is_decoy else entry_name
            self.psm_df.loc[i, 'gene'] = best_gene_name
            self.psm_df.loc[i, 'protein_description'] = best_description
            self.psm_df.loc[i, 'mapped_protein'] = ','.join(mapped_proteins)
            self.psm_df.loc[i, 'mapped_gene'] = ','.join(mapped_genes)
        (self.report_directory / 'peptide.pep.xml').unlink(missing_ok=True)
        (self.report_directory / 'interact.prot.xml').unlink(missing_ok=True)


    def add_app_score(self):
        app_score_paths = self.report_directory.glob('app_prediction.*.tsv')
        app_df = pd.DataFrame()
        app_df['sequence'] = np.unique(self.psm_df['sequence'])
        for app_score_path in app_score_paths:
            predictor = app_score_path.stem.split('.')[1]
            psm_app_df = pd.read_csv(app_score_path, sep='\t')
            if predictor == 'netmhcpan':
                seq_app_df = psm_app_df.loc[psm_app_df.groupby('Peptide')['EL_Rank'].idxmin(), ['Peptide', 'Allele', 'EL_Rank']]
                seq_app_df['netmhcpan_binder'] = seq_app_df['EL_Rank'].apply(lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            if predictor == 'mhcflurry':
                seq_app_df = psm_app_df.loc[psm_app_df.groupby('peptide')['mhcflurry_affinity_percentile'].idxmin(), ['peptide', 'allele', 'mhcflurry_affinity_percentile']]
                seq_app_df['mhcflurry_binder'] = seq_app_df['mhcflurry_affinity_percentile'].apply(lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            if predictor == 'bigmhc':
                seq_app_df = psm_app_df.loc[psm_app_df.groupby('pep')['BigMHC_EL'].idxmin(), ['pep', 'mhc', 'BigMHC_EL']]
                seq_app_df['BigMHC_EL'] *= 100
                seq_app_df['bigmhc_binder'] = seq_app_df['BigMHC_EL'].apply(lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            if predictor == 'netmhciipan':
                seq_app_df = psm_app_df.loc[psm_app_df.groupby('Peptide')['EL_Rank'].idxmin(), ['Peptide', 'Allele', 'EL_Rank']]
                seq_app_df['netmhciipan_binder'] = seq_app_df['EL_Rank'].apply(lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            if predictor == 'mixmhc2pred':
                psm_app_df['BestAllele'] = psm_app_df['BestAllele'].fillna('')
                psm_app_df['%Rank_best'] = psm_app_df['%Rank_best'].fillna(100)
                seq_app_df = psm_app_df.loc[psm_app_df.groupby('Peptide')['%Rank_best'].idxmin(), ['Peptide', 'BestAllele', '%Rank_best']]
                seq_app_df['mixmhc2pred_binder'] = seq_app_df['%Rank_best'].apply(lambda r: 'Strong' if r < 0.5 else ('Weak' if r < 2 else 'Non-binder'))
            seq_app_df.columns = ['sequence', 'best_allele', 'min_rank'] + [seq_app_df.columns[3]]
            app_df = pd.merge(app_df, seq_app_df, on='sequence', how='left', suffixes=('_left', '_right'))
            if 'min_rank_left' in app_df.columns:
                app_df['min_rank'] = app_df[['min_rank_left', 'min_rank_right']].min(axis=1)
                app_df['best_allele'] = np.where(app_df['min_rank'] == app_df['min_rank_left'], app_df['best_allele_left'], app_df['best_allele_right'])
                app_df.drop(columns=['min_rank_left', 'min_rank_right', 'best_allele_left', 'best_allele_right'], inplace=True)
        if 'best_allele' in app_df.columns:
            app_df['best_allele'] = app_df['best_allele'].apply(lambda a: normalize_allele_name(a.replace('__', '-')))
            app_df['binder'] = (app_df[[col for col in app_df.columns if '_binder' in col]]
                                .apply(lambda b: 'Strong' if 'Strong' in b.values else ('Weak' if 'Weak' in b.values else 'Non-binder'), axis=1))
        self.psm_df = pd.merge(self.psm_df, app_df, on='sequence', how='left')


    def generate_psm_report(self, psm_fdr=1, remove_decoy=False):
        if remove_decoy:
            psm_df = self.psm_df[self.psm_df['label'] == 'Target']
        else:
            psm_df = self.psm_df
        psm_df = psm_df[psm_df['psm_qvalue'] <= psm_fdr]
        psm_df.to_csv(self.report_directory / f'psm.tsv', index=False, sep='\t')
        return psm_df


    def generate_peptide_report(self, pep_fdr=1, remove_decoy=False):
        if remove_decoy:
            psm_df = self.psm_df[self.psm_df['label'] == 'Target']
        else:
            psm_df = self.psm_df
        psm_df = psm_df[psm_df['pep_qvalue'] <= pep_fdr]
        agg_dict = {
            'sequence': 'first',
            'prev_AA': lambda x: ','.join(set(x)),
            'next_AA': lambda x: ','.join(set(x)),
            'label': lambda x: 'Target' if 'Target' in set(x) else 'Decoy',
            'seq_len': 'first',
            'charge': lambda x: ','.join([str(c) for c in sorted(set(x))]),
            'score': 'max',
            'pep_qvalue': 'max',
            'protein': 'first',
            'protein_id': 'first',
            'entry_name': 'first',
            'protein_description': 'first',
            'mapped_protein': 'first',
            'mapped_gene': 'first'
        }
        for col in psm_df.columns:
            if 'binder' in col or 'best_allele' in col or 'min_rank' in col:
                agg_dict[col] = 'first'
        peptide_df = psm_df.groupby('peptide', as_index=False).agg(agg_dict)
        pep_stat = psm_df['peptide'].value_counts().reset_index().rename(columns={'count': 'spectral_count'})
        self.peptide_df = peptide_df.merge(pep_stat, how='left', on='peptide')
        self.peptide_df.to_csv(self.report_directory / f'peptide.tsv', index=False, sep='\t')
        return self.peptide_df


    def generate_sequence_report(self, seq_fdr=1, remove_decoy=False):
        if remove_decoy:
            psm_df = self.psm_df[self.psm_df['label'] == 'Target']
        else:
            psm_df = self.psm_df
        psm_df = psm_df[psm_df['seq_qvalue'] <= seq_fdr]
        agg_dict = {
            'peptide': lambda x: ','.join(set(x)),
            'prev_AA': lambda x: ','.join(set(x)),
            'next_AA': lambda x: ','.join(set(x)),
            'label': lambda x: 'Target' if 'Target' in set(x) else 'Decoy',
            'seq_len': 'first',
            'charge': lambda x: ','.join([str(c) for c in sorted(set(x))]),
            'score': 'max',
            'seq_qvalue': 'max',
            'protein': 'first',
            'protein_id': 'first',
            'entry_name': 'first',
            'protein_description': 'first',
            'mapped_protein': 'first',
            'mapped_gene': 'first'
        }
        for col in psm_df.columns:
            if 'binder' in col or 'best_allele' in col or 'min_rank' in col:
                agg_dict[col] = 'first'
        sequence_df = psm_df.groupby('sequence', as_index=False).agg(agg_dict)
        seq_stat = psm_df['sequence'].value_counts().reset_index().rename(
            columns={'count': 'spectral_count'})
        self.sequence_df = sequence_df.merge(seq_stat, how='left', on='sequence')
        self.sequence_df.to_csv(self.report_directory / f'sequence.tsv', index=False, sep='\t')
        return self.sequence_df


    def draw_result_figure(self, val_loss_list, psm_fdr=0.01):
        roc = calculate_roc(self.psm_df['psm_qvalue'], self.psm_df['label'])
        fig = plt.figure(constrained_layout=True, figsize=(10, 10))
        fig.suptitle(self.file_name, fontsize=16)
        gs = GridSpec(2, 2, figure=fig)
        colormap = get_cmap("tab10")

        final = fig.add_subplot(gs[0, 0])
        final.plot(*roc, c=colormap(0), ms='3', ls='none', marker='.', alpha=0.6)
        n_psms_at_fdr = np.sum((self.psm_df['psm_qvalue'] <= psm_fdr) & (self.psm_df['label'] == 'Target'))
        final.vlines(psm_fdr, 0, n_psms_at_fdr, ls='--', lw=1, color='k', alpha=0.7)
        final.hlines(n_psms_at_fdr, 0, psm_fdr, ls='--', lw=1, color='k', alpha=0.7)
        final.set_xlim((0, max(0.05, 2 * psm_fdr)))
        final.set_title('Final q-values')
        final.set_xlabel('q-value')
        final.set_ylabel('PSMs')
        final.set_ylim((0, final.get_ylim()[1]))

        dist = fig.add_subplot(gs[0, 1])
        scores = self.psm_df['score']
        labels = self.psm_df['label']
        _, bins, _ = dist.hist(scores[labels == 'Target'], label='Target', bins=30, alpha=0.5, color='g')
        dist.hist(scores[labels == 'Decoy'], label='Decoy', bins=bins, alpha=0.5, zorder=100, color='r')
        dist.set_title('Prediction distributions')
        dist.set_xlabel('Target probability')
        dist.set_ylabel('PSMs')
        dist.legend()

        loss = fig.add_subplot(gs[1, 1])
        min_x = []
        min_y = []
        for i, val_loss in enumerate(val_loss_list):
            loss.plot(range(1, len(val_loss) + 1), val_loss, c=colormap(i), marker=None, label=f'split {i+1}')
            min_y.append(np.min(val_loss))
            min_x.append(np.argmin(val_loss) + 1)
        loss.plot(min_x, min_y, ls='none', marker='x', ms='12', c='k', label='best models')
        loss.set_title('Validation loss')
        loss.set_xlabel('Epoch')
        loss.set_ylabel('Loss')
        loss.legend()
        plt.tight_layout()

        pdf_file = self.report_directory / 'training_report.pdf'
        pdf = plt_pdf.PdfPages(str(pdf_file), keep_empty=False)
        pdf.savefig(fig)
        pdf.close()
        plt.close(fig)

if __name__ == '__main__':
    psm_df = pd.read_csv('/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/Search_0225/JY_Class1_1M_DDA_60min_Slot1-10_1_541/psm.tsv', sep='\t')
    run_reporter = RunReporter(report_directory='/mnt/d/workspace/mhc-booster/experiment/JY_1_10_25M/Search_0225/JY_Class1_1M_DDA_60min_Slot1-10_1_541',
                               file_name='test', decoy_prefix='rev_')
    # psm_df['protein_id'] = ''
    # psm_df['entry_name'] = ''
    # psm_df['protein_description'] = ''
    # psm_df['mapped_protein'] = ''
    run_reporter.psm_df = psm_df
    # run_reporter.add_app_score()
    run_reporter.infer_protein('/mnt/d/data/Library/2025-02-25-decoys-contam-JY_var_splicing.fasta.fas')
    # run_reporter.generate_psm_report()
    # run_reporter.generate_peptide_report()
    # run_reporter.generate_sequence_report()
