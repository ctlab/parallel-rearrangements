from ete3 import SeqMotifFace
from collections import defaultdict, Counter
from itertools import chain

import numpy as np


def get_neighbour_motifs(n, ns_colors, offset, inverse=False):
    block, orient = n[:-1], n[-1]
    clr = ns_colors[block]

    if (orient == 'h' and not inverse) or (inverse and orient == 't'):
        return [[offset, offset + 30, "[]", None, 10, clr, clr, f"arial|2|white|+{block}"],
                [offset + 30, offset + 40, ">", None, 10, clr, clr, None]]
    else:
        return [[offset, offset + 10, "<", None, 10, clr, clr, None],
                [offset + 10, offset + 40, "[]", None, 10, clr, clr, f"arial|2|white|-{block}"]]


def generate_neighbour_face(node_ns, ns_colors, block, offsets):
    motifs = [[0, 0, "blank", None, 10, None, None, None]]

    for cur_offset, next_offset, node_n in zip(offsets, offsets.tolist()[1:], node_ns):
        if len(node_n) > 0:
            (nbr1, nbr2, copies, orientations) = node_n

            motifs.extend(get_neighbour_motifs(nbr1, ns_colors, cur_offset))

            for i, or_ in enumerate(orientations):
                motifs.extend(
                    get_neighbour_motifs(f'{block}' + ('h' if or_ == '+' else 't'), ns_colors, cur_offset + (i + 1) * 50))

            motifs.extend(get_neighbour_motifs(nbr2, ns_colors, cur_offset + (copies + 1) * 50, True))

            motifs.append([cur_offset + (copies + 1) * 50 + 40, next_offset, "blank", None, 10, None, None, None])
        else:
            motifs.append([cur_offset, next_offset, "blank", None, 10, None, None, None])

    for offset in offsets[1:-1]:
        motifs.append([offset - 10, offset - 10, "[]", None, 10, 'grey', 'grey', None])
        motifs.append([offset - 10, offset, "blank", None, 10, None, None, None])

    return SeqMotifFace('', motifs=motifs)


def align_neighbours(neighbours, all_genomes):
    aligned_neighbours = {}

    sorted_neighbours_keys = sorted(all_genomes, key=lambda k: len(neighbours[k]), reverse=True)
    max_blocks = len(neighbours[sorted_neighbours_keys[0]])

    l_neighbour_to_position = defaultdict(Counter)
    r_neighbour_to_position = defaultdict(Counter)

    align_smart = False
    for genome in sorted_neighbours_keys:
        node_ns = neighbours[genome]

        node_ns = sorted(node_ns, key=lambda k: (int(k[0][:-1]), int(k[1][:-1])))

        n = len(node_ns)
        # setting new indexes to better alignment
        if align_smart:
            new_node_ns = [()] * max_blocks
            used_nodes = [False for _ in range(n)]
            used_indexes = [False for _ in range(max_blocks)]

            for i, node_n in enumerate(node_ns):
                n1, n2, _1, _2 = node_n

                for new_ind, _freq in sorted(chain(l_neighbour_to_position[n1[:-1]].most_common(),
                                                   map(lambda v: (v[0], v[1] * 0.9),
                                                       r_neighbour_to_position[n2[:-1]].most_common())),
                                             key=lambda v: v[1], reverse=True):
                    if not used_indexes[new_ind]:
                        new_node_ns[new_ind] = node_n
                        used_nodes[i] = True
                        used_indexes[new_ind] = True
                        break

            free_old_index = [i for i, used in enumerate(used_nodes) if not used]
            free_new_index = [i for i, used in enumerate(used_indexes) if not used]

            assert len(free_old_index) <= len(free_new_index)

            for oi, ni in zip(free_old_index, free_new_index):
                new_node_ns[ni] = node_ns[oi]

            node_ns = new_node_ns

        # remember positions
        for i, node_n in enumerate(node_ns):
            if len(node_n) > 0:
                (n1, n2, _1, _2) = node_n
                l_neighbour_to_position[n1[:-1]][i] += 1
                r_neighbour_to_position[n2[:-1]][i] += 1

        align_smart = True
        aligned_neighbours[genome] = node_ns

    return aligned_neighbours


def get_offsets(neighbours):
    copies_2d = [[node_n[2] if len(node_n) > 0 else 0 for node_n in node_ns]
                 for node_ns in neighbours.values()]

    copies = np.max(copies_2d, axis=0)
    basic_offset = np.arange(0, 160 * len(copies) + 1, 160)

    shifted_copies = np.ones(len(copies) + 1, dtype='i4')
    shifted_copies[1:] = copies

    return basic_offset + np.cumsum((shifted_copies - 1) * 50)
