from collections import defaultdict

import numpy as np

import csv
import os


def get_characters_which_chr(blocks_df):
    def make_label(chrs):
        return 'chromosomes: ' + ', '.join(map(str, sorted(chrs)))

    characters = []

    which_chrs = defaultdict(lambda: defaultdict(list))
    block_lens = defaultdict(list)

    for (block, genome), df_block_genome in blocks_df.groupby(['block', 'species']):
        which_chrs[block][genome] = df_block_genome.chr.tolist()
        block_lens[block].extend(df_block_genome.chr_end - df_block_genome.chr_beg)

    for block in sorted(blocks_df.block.unique()):
        used_chroms = set()
        genome_colors = defaultdict(int)

        labels = ['no copies'] \
                 + sorted(list(set(make_label(chrs) for chrs in which_chrs[block].values())))

        for genome, chromosomes in which_chrs[block].items():
            label = make_label(chromosomes)
            used_chroms.update(chromosomes)

            color = labels.index(label)
            genome_colors[genome] = color

        multichromo = ';'.join(map(str, sorted(list(used_chroms))))
        mean_length = np.mean(block_lens[block])
        std_length = np.std(block_lens[block])

        characters.append((block, multichromo, mean_length, std_length, genome_colors, labels))

    return characters


def get_characters_stats_which_chr(characters, tree_holder):
    ans = []
    for block, multichromo, mean_length, std_length, genome_colors, labels in characters:
        tree_holder.count_innovations_fitch(genome_colors)

        score_rear, count_rear, count_all_rear = tree_holder.count_parallel_rearrangements(skip_grey=False)
        ans.append([block, multichromo, mean_length, std_length, score_rear, count_rear, count_all_rear, count_all_rear <= 1])

    return ans


def write_stats_csv_which_chr(stats, stats_file):
    rows = [['block', 'chroms', 'mean_length', 'std_length', 'parallel_rear_score', 'number_of_inconsistent_colors',
             'number_of_parallel_events', 'tree_consistent']] + stats
    with open(stats_file, 'w') as f:
        wtr = csv.writer(f)
        wtr.writerows(rows)


def write_characters_csv_which_chr(characters, folder):
    os.makedirs(folder, exist_ok=True)

    for i, (block, _1, _2, _3, genome_colors, labels) in enumerate(characters):
        rows = [['strain', 'character_state', 'character_state_annotation']] + \
               [[strain, color, labels[color]] for strain, color in genome_colors.items()]

        with open(folder + f'block_{block}.csv', 'w') as f:
            wtr = csv.writer(f)
            wtr.writerows(rows)


def write_trees_which_chr(characters, folder, show_branch_support, tree_holder, colors):
    os.makedirs(folder, exist_ok=True)

    for i, (block, _1, _2, _3, genome_colors, labels) in enumerate(characters):
        tree_holder.count_innovations_fitch(genome_colors)
        tree_holder.draw(folder + f'block_{block}.pdf', legend_labels=labels,
                         show_branch_support=show_branch_support, colors=colors)