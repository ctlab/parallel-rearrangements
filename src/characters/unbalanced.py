import csv
import os

from collections import defaultdict


def get_characters_stats_unbalanced(blocks, characters, tree_holder):
    ans = []
    for block, genome_colors in zip(blocks, characters):
        tree_holder.count_innovations_fitch(genome_colors)

        score_rear, count_rear, count_all_rear = tree_holder.count_parallel_rearrangements(skip_grey=False)
        ans.append([block, score_rear, count_rear, count_all_rear])

    return ans


def write_stats_csv_unbalanced(stats, cls, stats_file):
    rows = [['block', 'cluster', 'parallel_rear_score', 'parallel_rear_unique_innovation_count',
             'parallel_rear_all_innovations_count']] + \
           [stat[0:1] + [cl] + stat[1:] for stat, cl in zip(stats, cls)]
    with open(stats_file, 'w') as f:
        wtr = csv.writer(f)
        wtr.writerows(rows)


def call_unique_characters_one_cluster(cur_characters, cur_stats):
    count_blocks = defaultdict(list)
    for char, stats in zip(cur_characters, cur_stats):
        block = stats[0]
        count_blocks[frozenset(char.items())].append(block)
    return count_blocks


def call_unique_characters(cls_chars, cls_stats):
    return [call_unique_characters_one_cluster(c, s) for c, s in zip(cls_chars, cls_stats)]


def write_characters_csv_unbalanced(unique_chars_list, folder):
    fill_length = len(str(len(unique_chars_list)))
    for cl, unique_chars in enumerate(unique_chars_list):
        cl_folder = folder + f'cluster_{str(cl).zfill(fill_length)}/'
        os.makedirs(cl_folder, exist_ok=True)
        for unique_char, unique_count_blocks in unique_chars.items():
            # unfrozen = dict(unique_char)
            rows = [['strain', 'character_state', 'character_state_annotation']] + \
                   [[strain, color, f'{color} copies'] for strain, color in unique_char]

            with open(cl_folder + f'block{"s" if len(unique_count_blocks) > 1 else ""}'
                                  f'_{",".join(map(str,unique_count_blocks))}.csv', 'w') as f:
                wtr = csv.writer(f)
                wtr.writerows(rows)

def write_trees_unbalanced(unique_chars_list, folder, show_branch_support, tree_holder, colors):
    fill_length = len(str(len(unique_chars_list)))
    for cl, unique_chars in enumerate(unique_chars_list):
        cl_folder = folder + f'cluster_{str(cl).zfill(fill_length)}/'
        os.makedirs(cl_folder, exist_ok=True)

        for unique_char, unique_count_blocks in unique_chars.items():
            unfrozen = dict(unique_char)
            tree_holder.count_innovations_fitch(unfrozen)

            labels = [f'{i}{"+" if i == len(colors) - 1 else ""} copies'
                      for i in range(max(unfrozen.values()) + 1)]
            filepath = cl_folder + f'block{"s" if len(unique_count_blocks) > 1 else ""}_' \
                                   f'{",".join(map(str,unique_count_blocks))}.pdf'
            tree_holder.draw(filepath, legend_labels=labels, show_branch_support=show_branch_support, colors=colors)