from bg.grimm import GRIMMReader
from collections import Counter, defaultdict


class Unique_Filter:
    def __init__(self):
        self.blocks_copies = defaultdict(lambda: defaultdict(int))
        self.first_call = True

    def update_allowed_blocks(self, ps, strain):
        vs = [p[1] for p in ps]
        for v in vs:
            self.blocks_copies[v][strain] += 1

    def count_allowed(self, min_genomes):
        self.allowed_blocks = [b for b, cs in self.blocks_copies.items()
                               if all(map(lambda c: c == 1, cs.values())) and len(cs) >= min_genomes]

    def filter_unique(self, ps):
        return [p for p in ps if p[1] in self.allowed_blocks]

def grimm_filter_unique_gene(in_file, out_file, block_rate):
    lines = open(in_file).read().split('\n')
    strains = set()

    # make unique blocks list
    i = 0
    flt = Unique_Filter()
    while i < len(lines):
        line = lines[i]
        if GRIMMReader.is_genome_declaration_string(line):
            strain, chr = GRIMMReader.parse_genome_declaration_string(line).name.rsplit('.', 1)
            strains.add(strain)

            data_line = lines[i + 1]
            parsed = GRIMMReader.parse_data_string(data_line)[1]
            flt.update_allowed_blocks(parsed, strain)
            i += 2
        else:
            i += 1

    flt.count_allowed(block_rate * len(strains) / 100)
    # write allowed blocks
    i = 0
    with open(out_file, 'w') as f:
        while i < len(lines):
            line = lines[i]
            if GRIMMReader.is_genome_declaration_string(line):
                data_line = lines[i + 1]

                parsed = GRIMMReader.parse_data_string(data_line)[1]
                parsed = flt.filter_unique(parsed)

                print(line, file=f)
                print(' '.join(p[0] + p[1] for p in parsed), '@', file=f)
                i += 2
            else:
                i += 1

    return list(map(int, flt.allowed_blocks))

def filter_dataframe_unique(df):
    allowed_blocks = set()
    all_sp = len(df['species'].unique())
    for block, df_block in df.groupby('block'):
        if len(df_block) == len(df_block['species'].unique()) == all_sp:
            allowed_blocks.add(block)

    return df.loc[df['block'].isin(allowed_blocks)].copy()

def filter_dataframe_allowed(df, allowed_blocks):
    return df.loc[df['block'].isin(allowed_blocks)].copy()