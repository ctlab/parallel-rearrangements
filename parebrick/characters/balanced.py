from bg.grimm import GRIMMReader
from bg.tree import BGTree
from bg.genome import BGGenome
from bg.vertices import BGVertex


import numpy as np

import os
import csv

get_colors_by_edge = lambda e: e.multicolor.multicolors

def white_proportion(colors):
    return np.mean(list(map(lambda c: c == 0, colors)))

def get_character_by_edge(bg, edge, genomes, neighbour_index):
    def get_neighbour_with_genome(v, genome):
        return neighbour_index[(v, genome)]

    def get_genome_character_state_by_edge(genome):
        if cnt[BGGenome(genome)] == 1:
            return 0
        else:
            v1, v2 = edge.vertex1.name, edge.vertex2.name
            if v1 > v2: v1, v2 = v2, v1

            v1_neighbour = get_neighbour_with_genome(v1, genome)
            v2_neighbour = get_neighbour_with_genome(v2, genome)
            if bg.get_edge_by_two_vertices(v1_neighbour, v2_neighbour):
                pair = (v1_neighbour, v2_neighbour)
                if pair not in possible_edges:
                    possible_edges.append(pair)
                return 2 + possible_edges.index(pair)
            else:
                return 1

    cnt = get_colors_by_edge(edge)
    possible_edges = []
    return {genome: get_genome_character_state_by_edge(genome) for genome in genomes}, possible_edges

def construct_vertex_genome_index(bg):
    neighbour_index = {}
    for v in bg.bg:
        for edge in bg.get_edges_by_vertex(v):
            colors = get_colors_by_edge(edge)
            for color in colors:
                neighbour_index[(str(v), color.name)] = edge.vertex2

    return neighbour_index

def get_characters(grimm_file, genomes, logger):
    bg = GRIMMReader.get_breakpoint_graph(open(grimm_file))
    logger.info('Breakpoint graph parsed')

    logger.info(f'Edges in breakpoint graph: {len(list(bg.edges()))}')

    characters = []

    # consistency_checker = TreeConsistencyChecker(tree_file)
    for i, component_bg in enumerate(bg.connected_components_subgraphs()):
        nodes_len = len(list(component_bg.nodes()))
        if nodes_len == 2: continue

        logger.info(f'Getting characters from breakpoint graph component, size={len(component_bg.bg)}')

        neighbour_index = construct_vertex_genome_index(component_bg)

        for i_edge, edge in enumerate(component_bg.edges()):
            v1, v2 = edge.vertex1.name, edge.vertex2.name
            if v1 > v2: v1, v2 = v2, v1

            genome_colors, neighbour_edges = get_character_by_edge(component_bg, edge, genomes, neighbour_index)
            if white_proportion(genome_colors.values()) < 0.5: continue

            labels = ['edge exists', 'parallel edge doesn\'t exist'] + [f'inversion {v1}-{v2}-{v1n}-{v2n}'
                                                                        for (v1n, v2n) in neighbour_edges]

            characters.append((v1, v2, genome_colors, labels))

    return characters

def get_characters_stats_balanced(characters, tree_holder, distance_between_blocks):
    ans = []
    for v1, v2, genome_colors, labels in characters:
        tree_holder.count_innovations_fitch(genome_colors)

        b1, b2 = int(v1[:-1]), int(v2[:-1])
        if b1 > b2: b1, b2 = b2, b1

        score_rear, count_rear, count_all_rear = tree_holder.count_parallel_rearrangements(skip_grey=True)
        score_break, count_break = tree_holder.count_parallel_breakpoints()

        white_strains = [strain for strain, color in genome_colors.items() if color == 0]
        mean_break_length = np.mean([distance_between_blocks[(b1, b2)][strain] for strain in white_strains])

        ans.append([v1, v2, int(mean_break_length), score_rear, count_rear, count_all_rear, score_break, count_break])

    return ans

def write_stats_csv_balanced(stats, stats_file):
    rows = [['id', 'vertex1', 'vertex2', 'mean_break_length_nucleotide', 'parallel_rear_score', 'parallel_rear_unique_innovation_count',
            'parallel_rear_all_innovations_count', 'parallel_breakpoint_score', 'parallel_breakpoint_count']] + \
          [[i+1] + stat for i, stat in enumerate(stats)]
    with open(stats_file, 'w') as f:
        wtr = csv.writer(f)
        wtr.writerows(rows)

def write_characters_csv_balanced(characters, folder):
    os.makedirs(folder, exist_ok=True)
    fill_length = len(str(len(characters)))

    for i, (v1, v2, genome_colors, labels) in enumerate(characters):
        rows = [['strain', 'character_state', 'character_state_annotation']] + \
               [[strain, color, labels[color]] for strain, color in genome_colors.items()]

        with open(folder + f'id_{str(i + 1).zfill(fill_length)}_edge_{v1}-{v2}.csv', 'w') as f:
            wtr = csv.writer(f)
            wtr.writerows(rows)

def write_trees_balanced(characters, folder, show_branch_support, tree_holder, colors):
    os.makedirs(folder, exist_ok=True)
    fill_length = len(str(len(characters)))

    for i, (v1, v2, genome_colors, labels) in enumerate(characters):
        tree_holder.count_innovations_fitch(genome_colors)
        tree_holder.draw(folder + f'id_{str(i + 1).zfill(fill_length)}_edge_{v1}-{v2}.pdf', legend_labels=labels,
                         show_branch_support=show_branch_support, colors=colors)