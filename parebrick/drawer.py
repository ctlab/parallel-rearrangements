from parebrick.utils.data.parsers import parse_infercars_to_df
from parebrick.utils.data.stats import distance_between_blocks_distribution
from parebrick.utils.data.unique_gene_filters import filter_dataframe_unique

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import argparse
import os


def blocks_length_dist(df):
    return [end_ - start_ for start_, end_ in zip(df['chr_beg'], df['chr_end'])]


def lengths_between(state, log):
    plt.figure()
    if state == 'after':
        df_filtered = filter_dataframe_unique(df)
        ds = distance_between_blocks_distribution(df_filtered)
    else:
        ds = distance_between_blocks_distribution(df)

    sns.histplot(ds, bins=100, kde_kws={'log': log, 'alpha': 0.7})

    plt.ylabel('Number of blocks')
    plt.xlabel('Length in nucleotides')
    plt.title(f'Length of fragments not covered by {"common" if state == "after" else "any"} blocks')

    plt.tight_layout()
    plt.xlim(xmin=0)

    plt.savefig(out_folder + f'lengths_between_{state}_filtering{"_log" if log else ""}.pdf')
    # plt.show()


def number_of_genomes_weighted(weighted, log):
    plt.figure()
    vs, ws = [], []

    for _, df_block in df.groupby('block'):
        vs.append(len(df_block.species.unique()))
        lens = [row['chr_end'] - row['chr_beg'] for _, row in df_block.iterrows()]
        ws.append(np.mean(lens))

    bins = 100 if max(vs) > 100 else max(vs)
    if weighted:
        sns.histplot(vs, bins=bins, kde_kws={'log': log, 'alpha': 0.7, 'weights': ws})
    else:
        sns.histplot(vs, bins=bins, kde_kws={'log': log, 'alpha': 0.7})

    plt.ylabel('Length of fragments that are present\n in n genomes, nucleotides'
               if weighted else 'Number of blocks')
    plt.xlabel('Number of genomes')
    plt.title(f'{"Weighted f" if weighted else "F"}requency of synteny blocks')

    plt.xlim(xmin=0, xmax=max(vs))
    plt.tight_layout()
    plt.savefig(out_folder + f'blocks_frequency{"_weighted" if weighted else ""}{"_log" if log else ""}.pdf')
    # plt.show()


def block_length(log):
    plt.figure()
    df_filtered = filter_dataframe_unique(df)
    ds_before = blocks_length_dist(df)
    ds_after = blocks_length_dist(df_filtered)

    bins = np.linspace(min(ds_before), max(ds_before), 100)
    sns.histplot(ds_before, bins=bins, kde_kws={'log': log, 'alpha': 0.7}, label='not-common')
    sns.histplot(ds_after, bins=bins, kde_kws={'log': log, 'alpha': 0.7}, label='common', color='red')

    plt.ylabel('Number of blocks')
    plt.xlabel('Length in nucleotides')
    plt.xlim(xmin=0)
    plt.title(f'Distribution of synteny blocks length')

    plt.legend(loc=1)

    plt.tight_layout()
    plt.savefig(out_folder + f'block_lengths_distribution{"_log" if log else ""}.pdf')
    # plt.show()


def scatter_len_genomes_count(log):
    plt.figure()
    xs = []
    ys = []
    for block, df_block in df.groupby('block'):
        xs.append(len(df_block.species.unique()))
        lens = [row['chr_end'] - row['chr_beg'] for _, row in df_block.iterrows()]
        ys.append(np.mean(lens))

    if log: plt.yscale('log')
    plt.scatter(xs, ys, s=1)

    plt.xlabel('Number of genomes')
    plt.ylabel('Length of blocks')
    plt.title(f'Occurrence of synteny blocks vs its length')

    plt.xlim(xmin=0 - 0.02 * (max(xs)), xmax=max(xs) + 0.02 * (max(xs)))
    if not log: plt.ylim(ymin=0)

    plt.savefig(out_folder + f'scatter_number_length{"_log" if log else ""}.pdf')

    plt.subplots_adjust(left=0.14, right=0.99)

    plt.tight_layout()
    # plt.show()


def pan_blocks(permutations=10000):
    block_sets = [set(df_sp.block.unique()) for _, df_sp in df.groupby('species')]
    nbss, pnss = [], []
    for _ in range(permutations):
        block_sets = np.random.permutation(block_sets)
        nbs, pns = [], []
        accumulate_set = set()
        for bs in block_sets:
            left = bs - accumulate_set
            nbs.append(len(left))
            accumulate_set |= bs
            pns.append(len(accumulate_set))

        nbss.append(nbs)
        pnss.append(pns)

    nbss = np.array(nbss)
    xs = list(range(1, len(block_sets) + 1))

    def new_blocks():
        plt.figure()
        plt.plot(xs, np.median(nbss, axis=0), label='median (different permutations)')
        plt.fill_between(xs, np.percentile(nbss, 5, axis=0), np.percentile(nbss, 95, axis=0), alpha=0.4,
                         label='90% confidence interval')

        plt.xlabel('number of genomes')
        plt.ylabel('new blocks')
        plt.title('New blocks as a function of number of genomes')
        plt.legend(loc='upper right')

        plt.ylim(ymax=np.percentile(nbss, 95, axis=0)[1], ymin=0)
        plt.xlim(xmin=0, xmax=max(xs))

        # plt.subplots_adjust(top=0.99, right=0.99)
        plt.tight_layout()
        plt.savefig(out_folder + 'new_blocks.pdf')
        # plt.show()

    def pan():
        plt.figure()
        plt.plot(xs, np.median(pnss, axis=0), label='median (different permutations)', color='indianred')
        plt.fill_between(xs, np.percentile(pnss, 5, axis=0), np.percentile(pnss, 95, axis=0), alpha=0.4,
                         label='90% confidence interval', color='indianred')

        plt.xlabel('number of genomes')
        plt.ylabel('Pan-blocks count')
        plt.title('Pangenome')
        plt.legend(loc='lower right')

        plt.ylim(ymin=np.percentile(pnss, 5, axis=0)[0])
        plt.xlim(xmin=0, xmax=max(xs))
        # plt.subplots_adjust(top=0.99, right=0.99)
        plt.tight_layout()
        plt.savefig(out_folder + 'pan_blocks.pdf')
        # plt.show()

    new_blocks()
    pan()

def main():
    global out_folder, df

    parser = argparse.ArgumentParser(
        description='Building charts for pan-genome analysis based on synteny blocks.')

    parser.add_argument('--infercars_file', '-f', required=True,
                      help='Path to file in infercars format, can be found in main script output')

    parser.add_argument('--output', '-o', default='parebrick_charts', help='Path to output folder.')

    args = parser.parse_args()
    d = vars(args)

    file, out_folder = d['infercars_file'], d['output']

    if out_folder[-1] != '/': out_folder += '/'
    os.makedirs(out_folder, exist_ok=True)

    df = parse_infercars_to_df(file)
    sns.set(style="whitegrid", font="serif")

    print('Plotting lengths between blocks')
    lengths_between('before', log=False)
    lengths_between('after', log=False)
    lengths_between('before', log=True)
    lengths_between('after', log=True)

    print('Plotting number of genomes in blocks')
    number_of_genomes_weighted(weighted=False, log=False)
    number_of_genomes_weighted(weighted=True, log=False)
    number_of_genomes_weighted(weighted=False, log=True)
    number_of_genomes_weighted(weighted=True, log=True)

    print('Plotting blocks length distribution')
    block_length(log=True)
    block_length(log=False)

    print('Plotting scatter for occurrence of synteny blocks vs its length')
    scatter_len_genomes_count(False)
    scatter_len_genomes_count(True)

    print('Plotting pan-genome plots')
    pan_blocks()

if __name__ == "__main__":
    main()