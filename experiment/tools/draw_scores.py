
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def draw_scores(scores, labels, title = ''):
    bins = np.linspace(min(scores), max(scores), 100)
    plt.figure(figsize=(4, 3), dpi=200)
    plt.hist(scores[labels == 1], bins=bins, label='Target')
    plt.hist(scores[labels == -1], bins=bins, label='Decoy')
    plt.title(title)
    plt.xlabel('MS2 Entropy Score')
    plt.ylabel('Number of peptides')
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':


    pin_path = Path('/mnt/d/data/JY_1_10_25M/fragpipe/Search_best/JY_Class1_25M_DDA_60min_Slot1_12_1_552/JY_Class1_25M_DDA_60min_Slot1-12_1_552_edited.pin')
    pin_df = pd.read_csv(pin_path, sep='\t')
    draw_scores(pin_df['unweighted_spectral_entropy'].values, pin_df['Label'].values, 'JY in-house')
    pin_path = Path('/mnt/d/data/JPST002044/dda_pin/G230411_014_Slot2-28_1_9452_edited.pin')
    pin_df = pd.read_csv(pin_path, sep='\t')
    draw_scores(pin_df['unweighted_spectral_entropy'].values, pin_df['Label'].values, 'JPST002044')
    pin_path = Path('/mnt/d/data/PXD019643/fragpipe/Search_dda/160311_DK_AUT01-DN02_Bladder_W6-32_8_DDA_2_400-650mz_msms33_standard_edited.pin')
    pin_path = Path('/mnt/d/data/PXD019643/fragpipe/Search_dda/180122_AM_AUT01-DN06_Stomach_W6-32_10_DDA_3_400-650mz_msms27_standard_edited.pin')
    pin_df = pd.read_csv(pin_path, sep='\t')
    draw_scores(pin_df['unweighted_spectral_entropy'].values, pin_df['Label'].values, 'PXD019643 MSFragger DDA')
    pin_path = Path('/mnt/d/data/PXD019643/fragpipe/Search_ddaplus/160311_DK_AUT01-DN02_Bladder_W6-32_8_DDA_2_400-650mz_msms33_standard_edited.pin')
    pin_path = Path('/mnt/d/data/PXD019643/fragpipe/Search_ddaplus/180122_AM_AUT01-DN06_Stomach_W6-32_10_DDA_3_400-650mz_msms27_standard_edited.pin')
    pin_df = pd.read_csv(pin_path, sep='\t')
    draw_scores(pin_df['unweighted_spectral_entropy'].values, pin_df['Label'].values, 'PXD019643 MSFragger DDA+')
    # draw_scores(pin_df['delta_RT_loess'].values, pin_df['Label'].values)