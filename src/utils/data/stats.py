import numpy as np
import pandas as pd

from collections import defaultdict


def distance_between_blocks_distribution(df_blocks):
    ds = []
    for sp, df_sp in df_blocks.groupby('species'):
        df_sp = df_sp.sort_values(by=['chr_beg'])
        ds += (start_ - end_ for start_, end_ in zip(df_sp['chr_beg'][1:], df_sp['chr_end']))
    return ds


def distance_between_blocks_dict(df_blocks, genome_length):
    distances = defaultdict(lambda: defaultdict(lambda: 1e12))

    df_blocks['real_beg'] = df_blocks.apply(lambda row: row.chr_beg if row.orientation == '+' else row.chr_end, axis=1)
    df_blocks['real_end'] = df_blocks.apply(lambda row: row.chr_beg if row.orientation == '-' else row.chr_end, axis=1)

    for i, (strain, df_strain) in enumerate(df_blocks.groupby('species')):
        # print(strain, i, 'of', len(df_blocks['species'].unique()))
        l = genome_length[strain]

        xs = df_strain[['block', 'real_beg', 'real_end']].to_numpy()

        m, n = xs.shape

        combs = np.zeros((m, m, 2 * n), dtype=int)
        combs[:, :, :n] = xs[:, None, :]
        combs[:, :, n:] = xs
        combs.shape = (m * m, -1)

        combs = combs[np.where((combs[:, 3] - combs[:, 0]) > 0), :][0]

        for b1, st1, end1, b2, st2, end2 in combs:
            if b1 >= b2: continue
            distances[(b1, b2)][strain] = min(distances[(b1, b2)][strain],
                                                      abs(end1 - st2),
                                                      l - abs(end1 - st2),
                                                      abs(end2 - st1),
                                                      l - abs(end2 - st1)
                                                      )

    return distances


def check_stats_stains(tree, block_genomes):
    tree_genomes = tree.get_all_leafs()

    print('Strains in blocks file count:', len(block_genomes))
    print('Strains in tree leafs count:', len(tree_genomes))
    print('Intersected strains count:', len(block_genomes & tree_genomes))
    print()

    if not len(block_genomes & tree_genomes) == len(tree_genomes) == len(block_genomes):
        if len(block_genomes & tree_genomes) == 0:
            raise ValueError(
                'Seems like strains in block files and in tree leafs have different ids, intersection is empty.')
        left_tree = tree_genomes - (block_genomes & tree_genomes)
        if len(left_tree) > 0:
            print('PRUNING TREE')
            print('Those strains are thrown away from tree because there were not found in blocks data:',
                  ', '.join(left_tree))
            print()
            tree.prune(block_genomes & tree_genomes)
        left_blocks = block_genomes - (block_genomes & tree_genomes)
        if len(left_tree) > 0:
            print('FILTERING GENOMES IN BLOCKS')
            print('Those strains are thrown away from blocks because there were not found in tree leafs:',
                  ', '.join(left_blocks))
            print()

    return block_genomes & tree_genomes
