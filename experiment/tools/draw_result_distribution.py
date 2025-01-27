
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pyteomics import pepxml
from sympy.physics.vector.tests.test_printing import alpha


def get_charges(mhcv_path, netmhcpan_status_path):
    mhcv_df = pd.read_csv(mhcv_path, sep='\t')
    charges = []
    for peptide in mhcv_df['Peptide']:
        charges.append(peptide[-3])
    charges = np.array(charges, dtype=np.int32)

    netmhcpan_status_df = pd.read_csv(netmhcpan_status_path, sep='\t')
    non_binder_peptides = netmhcpan_status_df[netmhcpan_status_df['Binder'] == 'Non-binder']['Peptide']
    binder_peptides = netmhcpan_status_df[netmhcpan_status_df['Binder'] != 'Non-binder']['Peptide']
    non_binder_indices = mhcv_df['mhcv_peptide'].isin(non_binder_peptides)
    binder_indices = mhcv_df['mhcv_peptide'].isin(binder_peptides)

    return charges[binder_indices], charges[non_binder_indices]

def draw_charge_distribution(binder_charges, non_binder_charges):
    charge_labels = np.unique(np.concatenate((binder_charges, non_binder_charges)))
    binder_charge_count = [np.sum(binder_charges == i) for i in charge_labels]
    non_binder_charge_count = [np.sum(non_binder_charges == i) for i in charge_labels]
    plt.bar(charge_labels, binder_charge_count, width=0.4, color='blue', alpha=0.5, label='Binder')
    plt.bar(charge_labels + 0.4, non_binder_charge_count, width=0.4, color='orange', alpha=0.5, label='Non-binder')
    plt.xticks(charge_labels + 0.2, charge_labels.astype(int))
    plt.legend()
    plt.show()

def get_intensities(netmhcpan_status_path, peptide_path):
    netmhcpan_status_df = pd.read_csv(netmhcpan_status_path, sep='\t')
    non_binder_peptides = netmhcpan_status_df[netmhcpan_status_df['Binder'] == 'Non-binder']['Peptide']
    peptide_df = pd.read_csv(peptide_path, sep='\t')
    non_binder_indices = peptide_df['Peptide'].isin(non_binder_peptides)
    peptide_non_binder_df = peptide_df[non_binder_indices]
    peptide_binder_df = peptide_df[~non_binder_indices]
    return peptide_binder_df['Intensity'].to_numpy(dtype=np.float32), peptide_non_binder_df['Intensity'].to_numpy(dtype=np.float32)

def draw_intensity_distribution(binder_ints, non_binder_ints):
    binder_ints = np.log10(binder_ints + 1)
    non_binder_ints = np.log10(non_binder_ints + 1)
    ints = np.concatenate((binder_ints, non_binder_ints))
    bins = np.arange(np.floor(ints.min()), np.ceil(ints.max()), step=1).astype(int)
    plt.hist(binder_ints, bins=bins, color='blue', alpha=0.5, label='Binder')
    plt.hist(non_binder_ints, bins=bins, color='orange', alpha=0.5, label='Non-binder')
    plt.xlabel('Log10 Intensity')
    plt.ylabel('# Peptide')
    plt.legend()
    plt.show()
    plt.hist(binder_ints, bins=bins, color='blue', alpha=0.5, label='Binder', density=True)
    plt.hist(non_binder_ints, bins=bins, color='orange', alpha=0.5, label='Non-binder', density=True)
    plt.xlabel('Log10 Intensity')
    plt.ylabel('Density')
    plt.legend()
    plt.show()

def get_ims(mhcv_path, netmhcpan_status_path, pepxml_path):
    mhcv_df = pd.read_csv(mhcv_path, sep='\t')
    netmhcpan_status_df = pd.read_csv(netmhcpan_status_path, sep='\t')

    non_binder_peptides = netmhcpan_status_df[netmhcpan_status_df['Binder'] == 'Non-binder']['Peptide']
    binder_peptides = netmhcpan_status_df[netmhcpan_status_df['Binder'] != 'Non-binder']['Peptide']
    non_binder_indices = mhcv_df['mhcv_peptide'].isin(non_binder_peptides)
    binder_indices = mhcv_df['mhcv_peptide'].isin(binder_peptides)

    binder_scanids = mhcv_df['SpecId'][binder_indices].apply(lambda x: x.split('.')[-2]).values
    non_binder_scanids = mhcv_df['SpecId'][non_binder_indices].apply(lambda x: x.split('.')[-2]).values

    pepxml_data = pepxml.read(pepxml_path)
    pepxml_df = pd.DataFrame(pepxml_data)

    pepxml_binder_df = pepxml_df[pepxml_df['spectrumNativeID'].isin(binder_scanids)]
    binders = pd.DataFrame()
    binders['mz'] = pepxml_binder_df['precursor_neutral_mass'] / pepxml_binder_df['assumed_charge'] + 1.0072
    binders['rt'] = pepxml_binder_df['retention_time_sec'] / 60
    binders['im'] = pepxml_binder_df['ion_mobility']
    binders['charge'] = pepxml_binder_df['assumed_charge']

    pepxml_non_binder_df = pepxml_df[pepxml_df['spectrumNativeID'].isin(non_binder_scanids)]
    non_binders = pd.DataFrame()
    non_binders['mz'] = pepxml_non_binder_df['precursor_neutral_mass'] / pepxml_non_binder_df['assumed_charge'] + 1.0072
    non_binders['rt'] = pepxml_non_binder_df['retention_time_sec'] / 60
    non_binders['im'] = pepxml_non_binder_df['ion_mobility']
    non_binders['charge'] = pepxml_non_binder_df['assumed_charge']

    return binders, non_binders

def draw_ims_distribution(binders, non_binders):
    plt.figure(dpi=300)
    charges = np.unique(np.concatenate((binders['charge'], non_binders['charge'])))
    cmap = plt.get_cmap()
    colors = [cmap(i / (len(charges) - 1)) for i in range(len(charges))]
    binder_charges = np.unique(binders['charge'])
    for i, charge in enumerate(binder_charges):
        charge_indices = binders['charge'] == charge
        plt.scatter(binders['mz'][charge_indices].astype(float), binders['im'][charge_indices].astype(float), edgecolor='none', alpha=0.5, s=8, color=colors[i], label=charge)
    plt.legend()
    plt.show()

    plt.figure(dpi=300)
    non_binder_charges = np.unique(non_binders['charge'])
    for i, charge in enumerate(non_binder_charges):
        charge_indices = non_binders['charge'] == charge
        plt.scatter(non_binders['mz'][charge_indices].astype(float), non_binders['im'][charge_indices].astype(float), edgecolor='none', alpha=0.5, s=8, color=colors[i], label=charge)
    plt.legend()
    plt.show()

    plt.hist(binders['rt'], bins=20, color='blue', alpha=0.5, label='Binder', density=True)
    plt.hist(non_binders['rt'], bins=20, color='orange', alpha=0.5, label='Non-binder', density=True)
    plt.xlabel('RT')
    plt.ylabel('Density')
    plt.legend()
    plt.show()



if __name__ == '__main__':
    mhcv_path = '/mnt/d/workspace/mhc-validator-2/experiment/JY_1_10_25M/mhcvalidator_APP_PE_RT_MS2_CCS/JY_Class1_25M_DDA_60min_Slot1-12_1_552_edited_im2deep_MhcValidator/JY_Class1_25M_DDA_60min_Slot1-12_1_552_edited_im2deep.MhcValidator_annotated.tsv'
    netmhcpan_status_path = '/mnt/d/workspace/mhc-validator-2/experiment/JY_1_10_25M/mhcvalidator_APP_PE_RT_MS2_CCS/JY_Class1_25M_DDA_60min_Slot1-12_1_552_edited_im2deep_MhcValidator/JY_Class1_25M_DDA_60min_Slot1-12_1_552_edited_im2deep.netmhcpan_status.tsv'
    # binder_charges, non_binder_charges = get_charges(mhcv_path, netmhcpan_status_path)
    # draw_charge_distribution(binder_charges, non_binder_charges)
    #
    # peptide_path = '/mnt/d/data/JY_1_10_25M/fragpipe/Search_MSBooster_FDR1/JY_Class1_25M_DDA_60min_Slot1_12_1_552/peptide.tsv'
    # binder_ints, non_binder_ints = get_intensities(netmhcpan_status_path, peptide_path)
    # draw_intensity_distribution(binder_ints, non_binder_ints)

    pepxml_path = '/mnt/d/data/JY_1_10_25M/fragpipe/Search_MSBooster_FDR1/JY_Class1_25M_DDA_60min_Slot1_12_1_552/JY_Class1_25M_DDA_60min_Slot1-12_1_552.pepXML'
    binders, non_binders = get_ims(mhcv_path, netmhcpan_status_path, pepxml_path)
    draw_ims_distribution(binders, non_binders)
