import os
import re

import numpy as np
import pandas as pd

from itertools import takewhile
from io import StringIO
from bg.grimm import GRIMMReader
from collections import defaultdict, Counter

PATTERN = re.compile("([A-Za-z0-9_\(\)\/\s\.-]+)\.([A-Za-z0-9_]+):(\d+)-(\d+) ([+|-]).*")
COLUMNS = ["block", "species", "chr", "chr_beg", "chr_end", "orientation"]
BLOCKS_SEPARATOR = '-' * 80


def parse_infercars_to_df(file_name):
    def find_indices(lst, condition):
        return [i for i, elem in enumerate(lst) if condition(elem)]

    with open(file_name) as f:
        lines = f.readlines()

    last_line = len(lines) - 1
    while lines[last_line] == '\n': last_line -= 1

    n_at_end = len(lines) - 1 - last_line
    for _ in range(1 - n_at_end): lines.append('\n')

    bs = np.split(lines, find_indices(lines, lambda x: x[0] == ">"))
    temp = []

    for i, b in enumerate(bs):
        if len(b) == 0: continue
        b_i = int(b[0][1:])

        for oc in b[1:-1]:
            m = PATTERN.match(oc)
            temp.append([b_i, m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)])

    return pd.DataFrame(temp, columns=COLUMNS)


def export_df_to_infercars(df, file_name):
    with open(file_name, 'w') as f:
        for block, block_df in df.groupby('block'):
            print(f'>{block}', file=f)
            for i, row in block_df.iterrows():
                print(f'{row["species"]}.{row["chr"]}:{row["chr_beg"]}-{row["chr_end"]} {row["orientation"]}', file=f)
            print(file=f)


def genome_lengths_from_block_coords(in_file):
    with open(in_file) as f:
        head_lines = list(takewhile(lambda line: (line != BLOCKS_SEPARATOR + os.linesep) and
                                                 (line != BLOCKS_SEPARATOR + '\n'), f))

    # names of chromosomes
    df_head = pd.read_csv(StringIO(''.join(head_lines)), sep='\t')
    return {row['Description']: row['Size'] for index, row in df_head.iterrows()}


def get_genomes_contain_blocks_grimm(grimm_file):
    genomes, blocks = set(), set()

    with open(grimm_file) as f:
        ls = f.readlines()
    block_genome_count = defaultdict(Counter)

    for i in range(0, len(ls), 2):
        name = GRIMMReader.parse_genome_declaration_string(ls[i]).name
        data = GRIMMReader.parse_data_string(ls[i + 1])[1]
        genomes.add(name)
        for _, block in data:
            blocks.add(int(block))
            block_genome_count[int(block)][name] += 1

    return list(sorted(genomes)), list(sorted(blocks)), block_genome_count


def get_block_neighbours(grimm_file):
    block_neighbours = defaultdict(lambda: defaultdict(list))
    with open(grimm_file) as f: ls = f.readlines()

    i = 0
    while i < len(ls):
        l = ls[i]
        if GRIMMReader.is_genome_declaration_string(l):
            genome = GRIMMReader.parse_genome_declaration_string(l)
            data_line = ls[i + 1]
            bs = GRIMMReader.parse_data_string(data_line)[1]

            n = len(bs)
            j = 0

            while j < n:
                tandem_copies = 1
                prev_or, prev_block = bs[j % n]
                _, curr_block = bs[(j + 1) % n]
                next_or, next_block = bs[(j + 2) % n]

                if curr_block == prev_block:
                    j += 1
                    continue

                while curr_block == next_block:
                    j += 1
                    tandem_copies += 1
                    next_or, next_block = bs[(j + 2) % n]

                neighbours = (prev_block + ('h' if prev_or == '+' else 't'),
                              next_block + ('t' if next_or == '+' else 'h'))

                orientations = tuple(bs[(k + 1) % n][0] for k in range(j - tandem_copies + 1, j + 1))

                if orientations[0] == '-':
                    neighbours = (neighbours[1], neighbours[0])
                    orientations = tuple('+' if or_ == '-' else '+' for or_ in orientations[::-1])

                block_neighbours[int(curr_block)][genome.name].append((*neighbours, tandem_copies, orientations))

                j += 1

            i += 2
        else:
            i += 1

    return block_neighbours


def make_labels_dict(file, row_from='strain', row_to='label'):
    try:
        df = pd.read_csv(file)
        return {row[row_from]: row[row_to] for i, row in df.iterrows()}
    except FileNotFoundError:
        return None