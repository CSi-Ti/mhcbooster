import numpy as np
from typing import Union


def k_fold_split(s: np.array,
                 k_folds: int = 3,
                 random_state: Union[int, np.random.RandomState] = 0):

    if isinstance(random_state, int):
        random_state: np.random.RandomState = np.random.RandomState(random_state)
    else:
        random_state: np.random.RandomState = random_state

    indices = np.argsort(s)
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

    train_test_splits = list(zip(train_indices, val_indices))

    return train_test_splits