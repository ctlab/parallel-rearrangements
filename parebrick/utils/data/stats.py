import numpy as np

from collections import defaultdict


def distance_between_blocks_distribution(df_blocks):
    ds = []
    for sp, df_sp in df_blocks.groupby('species'):
        df_sp = df_sp.sort_values(by=['chr_beg'])
        ds += (start_ - end_ for start_, end_ in zip(df_sp['chr_beg'][1:], df_sp['chr_end']))
    return ds


def distance_between_blocks_dict(df_blocks, genome_length, allowed_blocks=None):
    distances = defaultdict(lambda: defaultdict(lambda: 1e12))

    if allowed_blocks is None:
        df = df_blocks.copy()
    else:
        df = df_blocks.loc[df_blocks['block'].isin(allowed_blocks)].copy()

    df['real_beg'] = df.apply(lambda row: row.chr_beg if row.orientation == '+' else row.chr_end, axis=1)
    df['real_end'] = df.apply(lambda row: row.chr_beg if row.orientation == '-' else row.chr_end, axis=1)

    for i, (strain, df_strain) in enumerate(df.groupby('species')):
        # print(strain, i, 'of', len(df['species'].unique()))
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
                                              abs(st1 - st2),
                                              l - abs(st1 - st2),
                                              abs(st1 - end2),
                                              l - abs(st1 - end2),
                                              abs(end1 - st2),
                                              l - abs(end1 - st2),
                                              abs(end1 - end2),
                                              l - abs(end1 - end2),
                                              )

    return distances


def check_stats_stains(tree, block_genomes, logger):
    tree_genomes = tree.get_all_leafs()

    logger.info(f'Strains in blocks file count: {len(block_genomes)}')
    logger.info(f'Strains in tree leafs count: {len(tree_genomes)}')
    logger.info(f'Intersected strains count: {len(block_genomes & tree_genomes)}')

    if not len(block_genomes & tree_genomes) == len(tree_genomes) == len(block_genomes):
        if len(block_genomes & tree_genomes) == 0:
            logger.error('Tree genomes: {tree_genomes}')
            logger.error('Blocks genomes: {block_genomes}')
            logger.error('Seems like strains in block files and in tree leafs have different ids, intersection is empty.')
            raise ValueError(
                'Seems like strains in block files and in tree leafs have different ids, intersection is empty.')
        left_tree = tree_genomes - (block_genomes & tree_genomes)
        if len(left_tree) > 0:
            logger.debug('PRUNING TREE')
            logger.debug(f'Those strains are thrown away from tree because there were not found in blocks data: {", ".join(left_tree)}')
            tree.prune(block_genomes & tree_genomes)
        left_blocks = block_genomes - (block_genomes & tree_genomes)
        if len(left_blocks) > 0:
            logger.debug('FILTERING GENOMES IN BLOCKS')
            logger.debug(f'Those strains are thrown away from blocks because there were not found in tree leafs: {", ".join(left_blocks)}')

    return block_genomes & tree_genomes


def get_coverages(df, genome_lengths):
    return [sum(df_strain['chr_end'] - df_strain['chr_beg']) / genome_lengths[strain]
            for strain, df_strain in df.groupby('species')]


def get_mean_coverage(df, genome_lengths):
    return np.mean(get_coverages(df, genome_lengths))
