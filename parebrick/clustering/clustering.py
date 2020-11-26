import numpy as np

from parebrick.clustering.distance_matrices import dist_matrix_similarity, dist_matrix_proximity

from sklearn.cluster import AgglomerativeClustering


def renew_index(cls):
    used, new_cls_indexes = set(), []
    for cl in cls:
        if not cl in used:
            used.add(cl)
            new_cls_indexes.append(cl)

    new_index_dict = {cl: i for i, cl in enumerate(new_cls_indexes)}
    return np.array([new_index_dict[cl] for cl in cls])

def clustering(characters, stats, distance_between_blocks, max_length, threshold, j, b, proximity_percentile, logger):
    J = dist_matrix_similarity(characters)
    if J.max() != 0: J /= J.max()
    logger.info('Jaccard index matrix constructed')

    B = dist_matrix_proximity(stats, distance_between_blocks, max_length, proximity_percentile)
    B /= B.max()
    logger.info('Proximity matrix constructed')

    D = J * j + B * b
    cls = AgglomerativeClustering(n_clusters=None, affinity='precomputed', linkage='average',
                                  distance_threshold=threshold).fit_predict(D)

    logger.info('Clustring is done')

    return renew_index(cls)


def split_by_cluster(characters, stats, cls):
    cls_chars, cls_stats = [], []
    for cl in np.unique(cls):
        inds = np.where(cls == cl)[0]
        cur_characters = [characters[i] for i in inds]
        cur_stats = [stats[i] for i in inds]

        cls_chars.append(cur_characters)
        cls_stats.append(cur_stats)
    return cls_chars, cls_stats