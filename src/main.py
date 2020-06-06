import argparse
import os

import numpy as np

from time import time
from itertools import takewhile
from operator import itemgetter

from src.characters.balanced import get_characters, write_characters_csv_balanced, get_characters_stats_balanced, \
    write_trees_balanced, write_stats_csv_balanced
from src.characters.unbalanced import get_characters_stats_unbalanced, write_stats_csv_unbalanced, \
    write_characters_csv_unbalanced, call_unique_characters, write_trees_unbalanced

from src.clustering.clustering import clustering, split_by_cluster

from src.utils.data.converters import block_coords_to_infercars
from src.utils.data.parsers import genome_lengths_from_block_coords, parse_infercars_to_df, \
    get_genomes_contain_blocks_grimm, make_labels_dict
from src.utils.data.unique_gene_filters import grimm_filter_unique_gene
from src.utils.data.stats import distance_between_blocks_dict, check_stats_stains

from src.utils.decorators import decorate

from src.tree.tree_holder import TreeHolder

# argument parsing
parser = argparse.ArgumentParser(
    description='Based on synteny blocks and phylogenetic tree this tool calls parallel rearrangements.')

parser.add_argument('--tree', '-t', required=True, help='Tree in newick format, must be parsable by ete3 library.'
                                                        'You can read more about formats supported by ete3 at http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees')
parser.add_argument('--blocks_folder', '-b', required=True,
                    help='Path to folder with blocks resulted as output of original Sibelia or maf2synteny tool.')
parser.add_argument('--labels', '-l', default='', help='Path to csv file with tree labels, must contain two columns: `strain` and `label`.')
parser.add_argument('--output', '-o', required=True, help='Path to output folder.')

show_branch_support = True
clustering_proximity_percentile = 25
clustering_threshold = 0.125
clustering_j = 0.75
clustering_b = 0.25

GRIMM_FILENAME = 'genomes_permutations.txt'
UNIQUE_GRIMM_FILENAME = 'genomes_permutations_unique.txt'

BLOCKS_COORD_FILENAME = 'blocks_coords.txt'
INFERCARS_FILENAME = 'blocks_coords.infercars'

STATS_FILE = 'stats.csv'
BALANCED_FOLDER = 'balanced_rearrangements_output/'
UNBALANCED_FOLDER = 'un' + BALANCED_FOLDER

CHARACTERS_FOLDER = 'characters/'
TREES_FOLDER = 'trees/'

BALANCED_COLORS = ['White', 'Gainsboro', 'LightGreen', 'LightBlue', 'NavajoWhite', 'LightPink', 'LightCoral', 'Purple',
                   'Navy', 'Olive', 'Teal', 'SaddleBrown', 'SeaGreen', 'DarkCyan', 'DarkOliveGreen', 'DarkSeaGreen']

UNBALANCED_COLORS = ['Gainsboro', 'White'] + BALANCED_COLORS[2:]

@decorate("Preprocess Data")
def preprocess_data():
    global unique_blocks
    unique_blocks = grimm_filter_unique_gene(blocks_folder + GRIMM_FILENAME, preprocessed_data_folder + UNIQUE_GRIMM_FILENAME)
    print('Converting block coords to infercars format')
    block_coords_to_infercars(blocks_folder + BLOCKS_COORD_FILENAME, preprocessed_data_folder + INFERCARS_FILENAME)


@decorate("Parsers and check strains")
def parsers_and_stats():
    global genome_lengths, blocks_df, tree_holder, genomes, blocks, block_genome_count

    genome_lengths = genome_lengths_from_block_coords(blocks_folder + BLOCKS_COORD_FILENAME)
    blocks_df = parse_infercars_to_df(preprocessed_data_folder + INFERCARS_FILENAME)
    tree_holder = TreeHolder(tree_file, labels_dict=make_labels_dict(labels_file))

    genomes, blocks, block_genome_count = get_genomes_contain_blocks_grimm(blocks_folder + GRIMM_FILENAME)

    genomes = check_stats_stains(tree_holder, set(genomes))


@decorate("Balanced rearrangements characters")
def balanced_rearrangements_characters():
    global b_characters
    b_characters = get_characters(preprocessed_data_folder + UNIQUE_GRIMM_FILENAME, genomes)

@decorate("Balanced rearrangements stats")
def balanced_rearrangements_stats():
    global b_characters, b_stats
    
    print('Counting distances between unique one-copy blocks, may take a while')
    distance_between_uniq_blocks = distance_between_blocks_dict(blocks_df, genome_lengths, unique_blocks)
    b_stats = get_characters_stats_balanced(b_characters, tree_holder, distance_between_uniq_blocks)
    char_stats = zip(b_characters, b_stats)

    print('Got characters after breakpoint graph consideration:', len(b_characters))
    # sorting for get most interesting characters in the beginning
    # [1][3] for vertex1,vertex2,parallel_rear_score, [1][6] for parallel_breakpoint_score, r[0][0] for vertex1
    char_stats = sorted(char_stats, key=lambda r: (r[1][3], r[1][6], r[0][0]), reverse=True)
    # remove convex characters
    char_stats = list(takewhile(lambda r: r[1][7] > 1, char_stats))

    # unzip
    b_characters = [char for char, stat in char_stats]
    b_stats = [stat for char, stat in char_stats]

    print('Left non-convex characters after filtering:', len(b_characters))

@decorate("Balanced rearrangements output")
def balanced_rearrangements_output():
    balanced_folder = output_folder + BALANCED_FOLDER
    os.makedirs(balanced_folder, exist_ok=True)

    stats_file = balanced_folder + STATS_FILE
    write_stats_csv_balanced(b_stats, stats_file)

    characters_folder = balanced_folder + CHARACTERS_FOLDER
    write_characters_csv_balanced(b_characters, characters_folder)

    trees_folder = balanced_folder + TREES_FOLDER
    write_trees_balanced(b_characters, trees_folder, show_branch_support, tree_holder, BALANCED_COLORS)

@decorate("Unbalanced rearrangements characters")
def unbalanced_rearrangements_characters():
    global ub_characters
    ub_characters = [{genome: block_genome_count[block][genome] for genome in genomes}
                     for block in blocks]

@decorate("Unbalanced rearrangements stats and clustering")
def unbalanced_rearrangements_stats_and_clustering():
    global ub_cls, ub_characters, ub_stats
    ub_stats = get_characters_stats_unbalanced(blocks, ub_characters, tree_holder)

    char_stats = zip(ub_characters, ub_stats)
    print('Got characters after copy number variation consideration:', len(blocks))
    # sorting for get most interesting characters in the beginning
    # [1][3] for vertex1,vertex2,parallel_rear_score, [1][6] for parallel_breakpoint_score, r[0][0] for vertex1
    char_stats = sorted(char_stats, key=lambda r: (r[1][1], r[1][0]), reverse=True)
    # remove convex characters
    char_stats = list(takewhile(lambda r: r[1][3] > 1, char_stats))

    print('Left non-convex characters after filtering:', len(char_stats))

    # unzip
    ub_characters = [char for char, stat in char_stats]
    ub_stats = [stat for char, stat in char_stats]

    print('Counting distances between non-convex character blocks, may take a while')
    distance_between_blocks = distance_between_blocks_dict(blocks_df, genome_lengths, set(map(itemgetter(0), ub_stats)))
    ub_cls = clustering(ub_characters, ub_stats, distance_between_blocks, max(genome_lengths.values()),
                        clustering_threshold, clustering_j, clustering_b, clustering_proximity_percentile)

    print('Clusters:', np.unique(ub_cls).shape[0])

@decorate('Unbalanced rearrangements output')
def unbalanced_rearrangements_output():
    unbalanced_folder = output_folder + UNBALANCED_FOLDER
    os.makedirs(unbalanced_folder, exist_ok=True)

    stats_file = unbalanced_folder + STATS_FILE
    write_stats_csv_unbalanced(ub_stats, ub_cls, stats_file)

    characters_folder = unbalanced_folder + CHARACTERS_FOLDER

    cls_chars, cls_stats = split_by_cluster(ub_characters, ub_stats, ub_cls)
    unique_chars_list = call_unique_characters(cls_chars, cls_stats)

    write_characters_csv_unbalanced(unique_chars_list, characters_folder)

    trees_folder = unbalanced_folder + TREES_FOLDER
    write_trees_unbalanced(unique_chars_list, trees_folder, show_branch_support, tree_holder, UNBALANCED_COLORS)


if __name__ == "__main__":
    start_time = time()
    d = vars(parser.parse_args())
    blocks_folder, output_folder, tree_file, labels_file = d['blocks_folder'], d['output'], d['tree'], d['labels']

    # folders
    if blocks_folder[:-1] != '/': blocks_folder += '/'
    if output_folder[:-1] != '/': output_folder += '/'

    preprocessed_data_folder = output_folder + 'preprocessed_data/'
    os.makedirs(preprocessed_data_folder, exist_ok=True)

    preprocess_data()
    parsers_and_stats()

    balanced_rearrangements_characters()
    balanced_rearrangements_stats()
    balanced_rearrangements_output()

    unbalanced_rearrangements_characters()
    unbalanced_rearrangements_stats_and_clustering()
    unbalanced_rearrangements_output()

    print('Total elapsed time:', time() - start_time)