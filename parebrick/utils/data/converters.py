import numpy as np
import pandas as pd

from io import StringIO

COLUMNS = ["block", "species", "chr", "chr_beg", "chr_end", "orientation"]
BLOCKS_SEPARATOR = '-' * 80


def block_coords_to_infercars(in_file, out_file):
    lines = open(in_file).read().split('\n')[:-1]
    ls = np.array(lines)
    bs = np.split(ls, np.where(ls == BLOCKS_SEPARATOR)[0])

    # names of chromosomes
    df_names = pd.read_csv(StringIO('\n'.join(bs[0])), sep='\t')
    chr_names = {}
    for index, row in df_names.iterrows():
        chr_names[row['Seq_id']] = row['Description']

    # blocks data
    with open(out_file, 'w') as f:
        for b in bs[1:-1]:
            block = b[1].split('#')[1]
            df_block = pd.read_csv(StringIO('\n'.join(b[2:])), sep='\t')

            print(f'>{block}', file=f)
            for index, row in df_block.iterrows():
                chr_name, start, end, strand = chr_names[row['Seq_id']], row['Start'], row['End'], row['Strand']
                print(chr_name + '.1:' + (f'{start}-{end}' if strand == '+' else f'{end}-{start}') + ' ' + strand, file=f)
            print(file=f)