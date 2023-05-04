import argparse
import os
import logging
import sys
import csv
import platform
import shutil

import numpy as np

from time import time
from itertools import takewhile
from operator import itemgetter

from parebrick.characters.balanced import get_characters_balanced, write_characters_csv_balanced, get_characters_stats_balanced, \
    write_trees_balanced, write_stats_csv_balanced
from parebrick.characters.unbalanced import get_characters_stats_unbalanced, write_stats_csv_unbalanced, \
    write_characters_csv_unbalanced, call_unique_characters, write_trees_unbalanced
from parebrick.characters.neighbours import write_trees_neightbours
from parebrick.characters.which_chromosome import get_characters_which_chr, get_characters_stats_which_chr, \
    write_stats_csv_which_chr, write_characters_csv_which_chr, write_trees_which_chr

from parebrick.clustering.clustering import clustering, split_by_cluster

from parebrick.utils.data.converters import block_coords_to_infercars
from parebrick.utils.data.parsers import genome_lengths_from_block_coords, parse_infercars_to_df, \
    get_genomes_contain_blocks_grimm, make_labels_dict, get_block_neighbours, export_df_to_infercars, \
    genome_genome_lengths_from_chromosomes_lengths
from parebrick.utils.data.unique_gene_filters import grimm_filter_unique_gene, filter_dataframe_unique, \
    filter_dataframe_allowed
from parebrick.utils.data.stats import distance_between_blocks_dict, check_stats_stains, get_mean_coverage

from parebrick.utils.decorators import decorate

from parebrick.tree.tree_holder import TreeHolder

logger = logging.getLogger()


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x


# argument parsing
def initialize():
    global parser, GRIMM_FILENAME, UNIQUE_GRIMM_FILENAME, BLOCKS_COORD_FILENAME, INFERCARS_FILENAME, STATS_FILE, \
        BALANCED_FOLDER, UNBALANCED_FOLDER, CHARACTERS_FOLDER, TREES_FOLDER, BALANCED_COLORS, UNBALANCED_COLORS, \
        clustering_proximity_percentile, clustering_threshold, clustering_j, clustering_j, clustering_b, \
        CSV_BLOCK_FILENAME, CSV_BLOCK_UNIQUE_FILENAME, CSV_GENOME_LENGTH, have_unique, NEIGHBOURS_FOLDER, \
        INFERCARS_UNIQUE_FILENAME, WHICH_CHR_FOLDER

    have_unique = True

    if platform.system() == 'Linux':
        os.environ['QT_QPA_PLATFORM'] = 'offscreen'

    parser = argparse.ArgumentParser(
        description='Based on synteny blocks and phylogenetic tree this tool calls parallel rearrangements '
                    '(not consistent with phylogenetic tree).')

    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')

    required.add_argument('--tree', '-t', required=True, help='Tree in newick format, must be parsable by ete3 library.'
                                                              'You can read more about formats supported by ete3 at http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees')
    required.add_argument('--blocks_folder', '-b', required=True,
                          help='Path to folder with blocks resulted as output of original Sibelia or maf2synteny tool.')

    optional.add_argument('--output', '-o', default='parebrick_output',
                          help='Path to output folder. Default: parebrick_output.')
    optional.add_argument('--labels', '-l', default='',
                          help='Path to csv file with tree labels, must contain two columns: `strain` and `label`.')

    optional.add_argument('--show_branch_support', '-sbs', type=str2bool, default=False, const=True, nargs='?',
                          help='Show branch support while tree rendering (ete3 parameter). Default: False.')

    optional.add_argument('--keep_non_parallel', '-knp', type=str2bool, default=True,
                          help='Keep rearrangements that are not parallel in result (consistent with phylogenetic tree). Default: True.')

    optional.add_argument('--filter_for_balanced', '-fb', type=float, default=80,
                          help='Minimal percentage of block occurrences in genomes for balanced rearrangements. '
                               'All blocks with lower occurrences rate will be removed. Default: 80.')

    optional.add_argument('--visualize_neighbours', '-vn', type=str2bool, default=True,
                          help='Use module for visualizing neighbours. Default: True.')


    optional.add_argument('--which_chr', '-chr', type=str2bool, default=False,
                          help='Use information about in which chromosome block is located. Default: False.')

    optional.add_argument('--clustering_tree_patterns_coef', '-j', type=restricted_float, metavar='[0-1]', default=0.8,
                          help='Coefficient (0-1) of weight of tree patterns similarity in clustering in unbalanced (copy number variation) module.'
                               'Rest of weight will be used in distance between blocks.'
                               'E.g. if set to 0.8 (default) distance between blocks coefficient will be set to 0.2.')

    optional.add_argument('--clustering_threshold', '-c', type=restricted_float, metavar='[0-1]', default=0.025,
                          help='Threshold for algorithm of clustering, default is 0.025.'
                               'Can be increased for getting larger clusters or decreased for getting smaller and more grouped clusters.')

    clustering_proximity_percentile = 25

    GRIMM_FILENAME = 'genomes_permutations.txt'
    UNIQUE_GRIMM_FILENAME = 'genomes_permutations_unique.txt'

    BLOCKS_COORD_FILENAME = 'blocks_coords.txt'
    INFERCARS_FILENAME = 'blocks_coords.infercars'
    INFERCARS_UNIQUE_FILENAME = 'blocks_unique_coords.infercars'

    CSV_BLOCK_FILENAME = 'blocks_coords.csv'
    CSV_BLOCK_UNIQUE_FILENAME = 'blocks_unique_coords.csv'
    CSV_GENOME_LENGTH = 'genomes_lengths.csv'

    STATS_FILE = 'stats.csv'
    BALANCED_FOLDER = 'balanced_rearrangements_output/'
    UNBALANCED_FOLDER = 'un' + BALANCED_FOLDER
    NEIGHBOURS_FOLDER = 'neighbours_output/'
    WHICH_CHR_FOLDER = 'which_chromosome_output/'

    CHARACTERS_FOLDER = 'characters/'
    TREES_FOLDER = 'tree_colorings/'

    BALANCED_COLORS = ['White', 'Gainsboro', 'DimGray', 'LightGreen', 'LightBlue', 'NavajoWhite', 'LightPink',
                       'DarkSeaGreen', 'Orchid', 'Navy', 'Olive', 'Teal', 'SaddleBrown', 'SeaGreen', 'DarkCyan',
                       'DarkOliveGreen', 'LightCoral']

    UNBALANCED_COLORS = ['Gainsboro', 'White'] + BALANCED_COLORS[3:]


# This module converts data from the output data format
# of Sibelia or Ragout scripts into the infercars format to simplify the subsequent annotation.
# Also filtering blocks in the grimm format for unique single-copy blocks for the breakpoint graph construction.
@decorate("Preprocess Data", logger)
def preprocess_data():
    global unique_blocks, balanced_block_rate
    unique_blocks = grimm_filter_unique_gene(blocks_folder + GRIMM_FILENAME,
                                             preprocessed_data_folder + UNIQUE_GRIMM_FILENAME, balanced_block_rate)
    logger.info('Converting block coords to infercars format')
    block_coords_to_infercars(blocks_folder + BLOCKS_COORD_FILENAME, preprocessed_data_folder + INFERCARS_FILENAME)


# In this module, parsing of input blocks and their coordinates, as well as a phylogenetic, takes place.
# Next, some statistics are calculated: the total number of blocks,
# the number of unique single-copy blocks, as well as the percentage of genome coverage with these two types of blocks.
# After that, the correspondence of strains set on the phylogenetic tree and strains set in the synthenic blocks
# is checked, if necessary, the missing strains are discarded.
@decorate("Parsers and check strains", logger)
def parsers_and_stats():
    global chr_lengths, blocks_df, tree_holder, genomes, blocks, block_genome_count, have_unique, unique_blocks, \
        neighbours

    chr_lengths = genome_lengths_from_block_coords(blocks_folder + BLOCKS_COORD_FILENAME)
    genome_lengths = genome_genome_lengths_from_chromosomes_lengths(chr_lengths)

    blocks_df = parse_infercars_to_df(preprocessed_data_folder + INFERCARS_FILENAME)
    unique_blocks_df = filter_dataframe_unique(blocks_df)
    filted_blocks_df = filter_dataframe_allowed(blocks_df, unique_blocks)
    neighbours = get_block_neighbours(blocks_folder + GRIMM_FILENAME)

    if len(filted_blocks_df) == 0:
        have_unique = False
        logger.warning('No unique one-copy blocks found. Balanced rearrangements will not be called')

    export_df_to_infercars(unique_blocks_df, preprocessed_data_folder + INFERCARS_UNIQUE_FILENAME)
    blocks_df.to_csv(preprocessed_data_folder + CSV_BLOCK_FILENAME, index=False)
    unique_blocks_df.to_csv(preprocessed_data_folder + CSV_BLOCK_UNIQUE_FILENAME, index=False)
    with open(preprocessed_data_folder + CSV_GENOME_LENGTH, 'w') as f:
        w = csv.writer(f)
        w.writerow(['Genome', 'Length'])
        w.writerows(chr_lengths.items())

    logger.info(f'Blocks count: {len(blocks_df["block"].unique())}')
    logger.info(f'Unique one-copy blocks count: {len(unique_blocks_df["block"].unique())}')

    logger.info(f'Mean block coverage: {get_mean_coverage(blocks_df, genome_lengths) * 100} %')
    logger.info(f'Mean common one-copy blocks coverage: '
                f'{get_mean_coverage(unique_blocks_df, genome_lengths) * 100 if have_unique else 0} %')
    logger.info(f'Mean {balanced_block_rate} % rate one-copy blocks coverage: '
                f'{get_mean_coverage(filted_blocks_df, genome_lengths) * 100 if have_unique else 0} %')

    tree_holder = TreeHolder(tree_file, logger, labels_dict=make_labels_dict(labels_file))

    genomes, blocks, block_genome_count = get_genomes_contain_blocks_grimm(blocks_folder + GRIMM_FILENAME)

    genomes = check_stats_stains(tree_holder, set(genomes), logger)


# In this module, the breakpoint of the graph is built using the bg library,
# then an algorithm for constructing characters and their states for each edge of the breakpoint graph is implemented.
@decorate("Balanced rearrangements characters", logger)
def balanced_rearrangements_characters():
    global b_characters
    b_characters = get_characters_balanced(preprocessed_data_folder + UNIQUE_GRIMM_FILENAME, genomes, logger)


# This module implements balanced rearrangements characters statistics calculation
# (measure for non-convexity proposed in the paper),
# as well as an estimate of the average distance of the length of the breakdown in the genome in nucleotides.
@decorate("Balanced rearrangements stats", logger)
def balanced_rearrangements_stats():
    global b_characters, b_stats, keep_consistent

    logger.info('Counting distances between unique one-copy blocks, may take a while')
    distance_between_uniq_blocks = distance_between_blocks_dict(blocks_df, chr_lengths, unique_blocks)
    b_stats = get_characters_stats_balanced(b_characters, tree_holder, distance_between_uniq_blocks)
    char_stats = zip(b_characters, b_stats)

    logger.info(f'Got characters after breakpoint graph consideration: {len(b_characters)}')
    # sorting for get most interesting characters in the beginning
    # [1][3] for vertex1,vertex2,parallel_rear_score, [1][6] for parallel_breakpoint_score, r[0][0] for vertex1
    char_stats = sorted(char_stats, key=lambda r: (r[1][2], r[1][5], r[0][0]), reverse=True)
    # remove convex characters
    if not keep_consistent:
        char_stats = list(takewhile(lambda r: not r[1][7], char_stats))
        logger.info(f'Left non-convex characters after filtering: {len(char_stats)}')
    else:
        logger.info('Percentage of consistent rearrangements (balanced): '
                    f'{np.mean([r[1][7] for r in char_stats]) * 100} %')

    # unzip
    b_characters = [char for char, stat in char_stats]
    b_stats = [stat for char, stat in char_stats]


# This module generates the output method files: the file stats.csv with the statistics
# of balanced rearrangements of all non-convex features.
# The tree_colorings folder contains characters visualized on the phylogenetic tree in .pdf format,
# different colors correspond to different states of the character, and the characters folder contains the .csv format
# files for the character states in different strains for each character.
@decorate("Balanced rearrangements output", logger)
def balanced_rearrangements_output():
    balanced_folder = output_folder + BALANCED_FOLDER
    os.makedirs(balanced_folder, exist_ok=True)

    stats_file = balanced_folder + STATS_FILE
    write_stats_csv_balanced(b_stats, stats_file)

    characters_folder = balanced_folder + CHARACTERS_FOLDER
    write_characters_csv_balanced(b_characters, characters_folder)

    trees_folder = balanced_folder + TREES_FOLDER
    write_trees_balanced(b_characters, trees_folder, show_branch_support, tree_holder, BALANCED_COLORS)


# In this module, the mapping of blocks by its number for each strain is performed.
@decorate("Unbalanced rearrangements characters", logger)
def unbalanced_rearrangements_characters():
    global ub_characters
    ub_characters = [{genome: block_genome_count[block][genome] for genome in genomes}
                     for block in blocks]


# This module implements balanced rearrangements characters statistics calculation
# (measure for non-convexity proposed in the paper),
# as well as the construction of distance matrices based on a measure of Jacquard similarity and distance in nucleotides
# and the subsequent clustering for unbalanced characters.
@decorate("Unbalanced rearrangements stats and clustering", logger)
def unbalanced_rearrangements_stats_and_clustering():
    global ub_cls, ub_characters, ub_stats
    ub_stats = get_characters_stats_unbalanced(blocks, ub_characters, tree_holder)

    char_stats = zip(ub_characters, ub_stats)
    logger.info(f'Got characters after copy number variation consideration: {len(blocks)}')
    # sorting for get most interesting characters in the beginning
    # [1][3] for vertex1,vertex2,parallel_rear_score, [1][6] for parallel_breakpoint_score, r[0][0] for vertex1
    char_stats = sorted(char_stats, key=lambda r: (r[1][1], r[1][0]), reverse=True)
    # remove convex characters
    if not keep_consistent:
        char_stats = list(takewhile(lambda r: not r[1][4], char_stats))
        logger.info(f'Left non-convex characters after filtering: {len(char_stats)}')
    else:
        logger.info('Percentage of consistent rearrangements (unbalanced): '
                    f'{np.mean([r[1][4] for r in char_stats]) * 100} %')

    # unzip
    ub_characters = [char for char, stat in char_stats]
    ub_stats = [stat for char, stat in char_stats]

    if len(ub_stats) == 0:
        logger.info('There is no non-convex characters, skipping clustering part')
        return
    logger.info('Counting distances between non-convex character blocks, may take a while')
    if len(ub_stats) > 1:
        distance_between_blocks = distance_between_blocks_dict(blocks_df, chr_lengths,
                                                               set(map(itemgetter(0), ub_stats)))
        ub_cls = clustering(ub_characters, ub_stats, distance_between_blocks, max(chr_lengths.values()),
                            clustering_threshold, clustering_j, clustering_b, clustering_proximity_percentile, logger)
    else:
        ub_cls = np.array([0])

    logger.info(f'Clusters: {np.unique(ub_cls).shape[0]}')


# In this module, the output of the method files is generated: the file stats.csv with the statistics of unbalanced
# genomic rearrangements of all non-convex characters.
# The tree_colorings folder contains characters visualized on the phylogenetic tree in .pdf format,
# and in the characters folder contain the .csv files for the character states in different strains
# for all the characters.
# These .pdf and .csv files are located in subfolders according to their presentation in clustering.
@decorate('Unbalanced rearrangements output', logger)
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


@decorate('Visualize neighbours output', logger)
def neighbours_output():
    neighbours_folder = output_folder + NEIGHBOURS_FOLDER
    os.makedirs(neighbours_folder, exist_ok=True)

    write_trees_neightbours([s[0] for s in ub_stats], ub_characters, neighbours, neighbours_folder, show_branch_support,
                            tree_holder, UNBALANCED_COLORS)


@decorate("Which chromosome characters", logger)
def which_chromosome_characters():
    global chr_characters
    chr_characters = get_characters_which_chr(blocks_df)


@decorate("Which chromosome stats and clustering", logger)
def which_chromosome_stats_and_clustering():
    global chr_cls, chr_characters, chr_stats
    chr_stats = get_characters_stats_which_chr(chr_characters, tree_holder)

    char_stats = zip(chr_characters, chr_stats)
    logger.info(f'Got characters after which chromosome consideration: {len(blocks)}')
    # sorting for get most interesting characters in the beginning
    char_stats = sorted(char_stats, key=lambda r: (r[1][4], r[1][3]), reverse=True)
    # remove convex characters
    if not keep_consistent:
        char_stats = list(takewhile(lambda r: not r[1][7], char_stats))
        logger.info(f'Left non-convex characters after filtering: {len(char_stats)}')
    else:
        logger.info('Percentage of consistent rearrangements (unbalanced): '
                    f'{np.mean([r[1][7] for r in char_stats]) * 100} %')

    # unzip
    chr_characters = [char for char, stat in char_stats]
    chr_stats = [stat for char, stat in char_stats]


@decorate('Which chromosome output', logger)
def which_chromosome_output():
    chr_folder = output_folder + WHICH_CHR_FOLDER
    os.makedirs(chr_folder, exist_ok=True)

    stats_file = chr_folder + STATS_FILE
    write_stats_csv_which_chr(chr_stats, stats_file)

    characters_folder = chr_folder + CHARACTERS_FOLDER

    write_characters_csv_which_chr(chr_characters, characters_folder)

    trees_folder = chr_folder + TREES_FOLDER
    write_trees_which_chr(chr_characters, trees_folder, show_branch_support, tree_holder, UNBALANCED_COLORS)


def main():
    def logger_initialize():
        file_handler = logging.FileHandler(filename=output_folder + 'run.log')
        stdout_handler = logging.StreamHandler(sys.stdout)
        handlers = [file_handler, stdout_handler]

        logging.basicConfig(
            level=logging.DEBUG,
            format='[%(asctime)s] %(levelname)s - %(message)s',
            handlers=handlers
        )

    global blocks_folder, output_folder, tree_file, labels_file, preprocessed_data_folder, show_branch_support, \
        have_unique, keep_consistent, balanced_block_rate, clustering_threshold, clustering_j, clustering_b
    initialize()

    start_time = time()
    d = vars(parser.parse_args())
    blocks_folder, output_folder, tree_file, labels_file, show_branch_support, keep_consistent, balanced_block_rate, \
    visualize_neighbours, clustering_j, clustering_threshold, which_chr_flag = d['blocks_folder'], d['output'], \
           d['tree'], d['labels'], d['show_branch_support'], d['keep_non_parallel'], d['filter_for_balanced'], \
           d['visualize_neighbours'], d['clustering_tree_patterns_coef'], d['clustering_threshold'], d['which_chr']

    clustering_b = 1 - clustering_j

    # folders
    if blocks_folder[-1] != '/': blocks_folder += '/'
    if output_folder[-1] != '/': output_folder += '/'
    shutil.rmtree(output_folder, ignore_errors=True)
    print('Cleaning output folder')

    preprocessed_data_folder = output_folder + 'preprocessed_data/'
    os.makedirs(preprocessed_data_folder, exist_ok=True)

    logger_initialize()

    preprocess_data()
    parsers_and_stats()

    if have_unique:
        balanced_rearrangements_characters()
        balanced_rearrangements_stats()
        balanced_rearrangements_output()

    unbalanced_rearrangements_characters()
    unbalanced_rearrangements_stats_and_clustering()
    if len(ub_stats) != 0: unbalanced_rearrangements_output()

    if which_chr_flag:
        which_chromosome_characters()
        which_chromosome_stats_and_clustering()
        which_chromosome_output()

    if visualize_neighbours:
        neighbours_output()

    logger.info(f'Total elapsed time: {time() - start_time} seconds')


if __name__ == "__main__":
    main()
