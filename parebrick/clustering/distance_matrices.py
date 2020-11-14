import numpy as np


def dist_matrix_similarity(characters):
    n = len(characters)
    arr = np.zeros((n, n))
    for i1, char1 in enumerate(characters):
        for i2, char2 in enumerate(characters):
            if i1 >= i2: continue
            arr[i1, i2] = arr[i2, i1] = \
                sum(char1[genome] != char2[genome] for genome in char1.keys() | char2.keys())

    return arr


def dist_matrix_proximity(stats, distance_between_blocks, max_length, percentile):
    n = len(stats)
    arr = np.zeros((n, n))
    for i1, stats1 in enumerate(stats):
        for i2, stats2 in enumerate(stats):
            block1, block2 = stats1[0], stats2[0]
            if i1 >= i2: continue
            b1, b2 = block1, block2
            if b1 > b2: b1, b2 = b2, b1
            ds = distance_between_blocks[(b1, b2)].values()
            arr[i1, i2] = arr[i2, i1] = np.percentile(list(ds), percentile) if len(ds) > 0 else max_length

    return arr
