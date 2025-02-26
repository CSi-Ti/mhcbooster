import pandas as pd

from src.predictors.netmhcpan_helper import NetMHCpanHelper
from src.predictors.mhcflurry_helper import MhcFlurryHelper
from src.predictors.bigmhc_helper import BigMhcHelper
from src.predictors.mixmhc2pred_helper import MixMhc2PredHelper

class BindingPredictor:
    def __init__(self, result_folder, allele_path, app_predictors):
        self.result_folder = result_folder
        self.app_predictors = [predictor.lower() for predictor in app_predictors]

    def predict(self):
        bind_df = pd.DataFrame()
        for run_folder in self.result_folder.iterdir():

            result_file = run_folder.glob('*.result.tsv') # TODO
            if len(result_file) == 0:
                print(f'Skipping... No result file found in {run_folder}')
                continue

            if 'netmhcpan' in self.app_predictors:
                netmhcpan_file = run_folder.glob('*.netmhcpan.tsv')
                if len(netmhcpan_file) == 0:

                    self.run_netmhcpan

                else:
                    netmhcpan_file = netmhcpan_file[0]


def run_netmhcpan(peptides, alleles, strong_rank=0.5, weak_rank=2):
    netmhcpan_helper = NetMHCpanHelper(peptides=peptides, alleles=alleles, mhc_class='I', n_threads=8)
    netmhcpan_helper.make_predictions()

    binding_status = []
    for pep in peptides:
        pep_pred = netmhcpan_helper.predictions[pep]
        highest_rank = 100
        for allele in netmhcpan_helper.alleles:
            allele_rank = pep_pred[allele]['el_rank']
            print(pep_pred[allele]['binder'])
            if allele_rank < highest_rank:
                highest_rank = allele_rank
        if highest_rank < strong_rank:
            binding_status.append('Strong')
        elif highest_rank < weak_rank:
            binding_status.append('Weak')
        else:
            binding_status.append('Non-binder')
    return binding_status


def run_mhcflurry(peptides, alleles, strong_rank=0.5, weak_rank=2):
    mhcflurry_helper = MhcFlurryHelper(peptides=peptides, alleles=alleles, report_directory=None)
    mhcflurry_helper.predict_df()

    binding_status = []
    for pep in peptides:
        pep_df = mhcflurry_helper.pred_df[mhcflurry_helper.pred_df['peptide'] == pep]
        rank = pep_df['mhcflurry_affinity_percentile'].min()
        if rank < strong_rank:
            binding_status.append('Strong')
        elif rank < weak_rank:
            binding_status.append('Weak')
        else:
            binding_status.append('Non-binder')
    return binding_status


def run_bigmhc(peptides, alleles, strong_rank=0.5, weak_rank=2):
    bigmhc_helper = BigMhcHelper(peptides=peptides, alleles=alleles)
    bigmhc_helper.predict_df()
    binding_status = []
    for pep in peptides:
        pep_df = bigmhc_helper.pred_df[bigmhc_helper.pred_df['pep'] == pep]
        rank = pep_df['BigMHC_EL'].min()
        if rank < strong_rank / 100.0:
            binding_status.append('Strong')
        elif rank < weak_rank / 100.0:
            binding_status.append('Weak')
        else:
            binding_status.append('Non-binder')
    return binding_status


def run_netmhciipan(peptides, alleles, strong_rank=0.5, weak_rank=2):
    netmhcpan_helper = NetMHCpanHelper(peptides=peptides, alleles=alleles, mhc_class='II', n_threads=8)
    netmhcpan_helper.make_predictions()

    binding_status = []
    for pep in peptides:
        pep_pred = netmhcpan_helper.predictions[pep]
        highest_rank = 100
        for allele in netmhcpan_helper.alleles:
            allele_rank = pep_pred[allele]['el_rank']
            print(pep_pred[allele]['binder'])
            print(allele_rank)
            if allele_rank < highest_rank:
                highest_rank = allele_rank
        if highest_rank < strong_rank:
            binding_status.append('Strong')
        elif highest_rank < weak_rank:
            binding_status.append('Weak')
        else:
            binding_status.append('Non-binder')
    return binding_status


def run_mixmhc2pred(peptides, alleles, strong_rank=0.5, weak_rank=2):
    mixmhc2pred_helper = MixMhc2PredHelper(peptides=peptides, alleles=alleles)
    mixmhc2pred_helper.predict_df()
    binding_status = []
    for pep in peptides:
        pep_df = mixmhc2pred_helper.pred_df[mixmhc2pred_helper.pred_df['Peptide'] == pep]
        rank = pep_df['%Rank_best'].values[0]
        if rank == 'NA':
            binding_status.append('Non-binder')
        else:
            rank = float(rank)
            if rank < strong_rank:
                binding_status.append('Strong')
            elif rank < weak_rank:
                binding_status.append('Weak')
            else:
                binding_status.append('Non-binder')
    return binding_status


def run_predictor(predictor, peptides, alleles, strong_rank=0.5, weak_rank=2):
    binding_status = ''
    if predictor == 'netmhcpan':
        binding_status = run_netmhcpan(peptides, alleles, strong_rank, weak_rank)
    if predictor == 'mhcflurry':
        binding_status = run_mhcflurry(peptides, alleles, strong_rank, weak_rank)
    if predictor == 'bigmhc':
        binding_status = run_bigmhc(peptides, alleles, strong_rank, weak_rank)
    if predictor == 'netmhciipan':
        binding_status = run_netmhciipan(peptides, alleles, strong_rank, weak_rank)
    if predictor == 'mixmhc2pred':
        binding_status = run_mixmhc2pred(peptides, alleles, strong_rank, weak_rank)
    return binding_status


if __name__ == '__main__':
    peptides = ['AAAAAAAAL', 'AAAAGRIAI', 'AAAAGRKTL', 'AAAAPALAAAA', 'AAAAPRAVA']
    # alleles = ['HLA-A0201', 'HLA-B0702', 'HLA-C0702']
    alleles = ['DRB1*04:01', 'DRB4*01:01']
    status = run_mixmhc2pred(peptides, alleles)
    print(status)