import numpy as np
from typing import Iterable, Union
from collections import Counter
from itertools import cycle


def k_fold_split(s: np.array,
                 peptides: Iterable[str],
                 k_folds: int = 3,
                 random_state: Union[int, np.random.RandomState] = 0):

    if isinstance(random_state, int):
        random_state: np.random.RandomState = np.random.RandomState(random_state)
    else:
        random_state: np.random.RandomState = random_state
    # get the counts of each peptide sequence
    peptide_counts = Counter(peptides)

    unique_peptide_indices = []
    nonunique_peptide_indices = {}

    # make lists of unique and non-unique peptide indices.
    for i, pep in enumerate(peptides):
        if peptide_counts[i] == 1:
            unique_peptide_indices.append(i)
        else:
            if pep not in nonunique_peptide_indices:
                nonunique_peptide_indices[pep] = []
            nonunique_peptide_indices[pep].append(i)

    unique_peptide_indices = np.array(unique_peptide_indices).astype(int)
    indices = unique_peptide_indices[np.argsort(s[unique_peptide_indices])]
    fold_indices = [[] for _ in range(k_folds)]
    k = list(range(k_folds))
    for i in range(len(indices) // k_folds):
        random_state.shuffle(k)
        for j in range(k_folds):
            fold_indices[k[j]].append(indices[i * k_folds + j])
    for i in range(len(indices) % k_folds):
        fold_indices[i].append(indices[len(indices) - 1 - i])

    train_indices = [[] for _ in range(k_folds)]
    val_indices = [[] for _ in range(k_folds)]
    for i in range(k_folds):
        for j in range(k_folds):
            if j != i:
                train_indices[i] += fold_indices[j]
        val_indices[i] = fold_indices[i]

    # randomly get one of the splits to start
    k = list(range(k_folds))
    random_state.shuffle(k)
    k = cycle(k)
    for indices in nonunique_peptide_indices.values():
        val_k = next(k)
        val_indices[val_k] += indices  # add the indices to a single validation set

        # and now add them to the training sets in all the other splits
        train_ks = list(range(k_folds))
        train_ks.remove(val_k)
        for i in train_ks:
            train_indices[i] += indices

    for k in range(k_folds):
        random_state.shuffle(train_indices[k])
        random_state.shuffle(val_indices[k])
        train_indices[k] = np.array(train_indices[k])
        val_indices[k] = np.array(val_indices[k])

    train_test_splits = list(zip(train_indices, val_indices))

    return train_test_splits