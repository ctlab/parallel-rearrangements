from bg.grimm import GRIMMReader
from collections import Counter


class Unique_Filter:
    def __init__(self):
        self.allowed_blocks = {}
        self.first_call = True

    def update_allowed_blocks(self, ps):
        vs = [p[1] for p in ps]
        if self.first_call == True:
            self.allowed_blocks = vs
            self.first_call = False
        counter = Counter(vs)
        self.allowed_blocks = {b for b in self.allowed_blocks if counter[b] == 1}

    def filter_unique(self, ps):
        return [p for p in ps if p[1] in self.allowed_blocks]

def grimm_filter_unique_gene(in_file, out_file):
    lines = open(in_file).read().split('\n')

    # make unique blocks list
    i = 0
    flt = Unique_Filter()
    while i < len(lines):
        line = lines[i]
        if GRIMMReader.is_genome_declaration_string(line):
            data_line = lines[i + 1]
            parsed = GRIMMReader.parse_data_string(data_line)[1]
            flt.update_allowed_blocks(parsed)
            i += 2
        else:
            i += 1

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

