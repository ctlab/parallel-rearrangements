# PaReBrick: PArallel REarrangements and BReakpoints identification toolkit

---
## Motivation
High plasticity of bacterial genomes is provided by numerous mechanisms including horizontal gene transfer and recombination via numerous flanking repeats. 
Genome rearrangements such as inversions, deletions, insertions, and duplications may independently occur in different strains, providing parallel adaptation. 
Specifically, such rearrangements might be responsible for multi-virulence, antibiotic resistance, and antigenic variation. 
However, identification of such events requires laborious manual inspection and verification of phyletic pattern consistency.

## Methods and Results

![Pipeline of tool](figs/pipeline.svg)

We present **tool `PaReBrick`** — implementation of an algorithmic solution for the identification of parallel rearrangements in bacterial population.
We define the term "parallel rearrangements" as events that occur independently in phylogenetically distant bacterial strains and present a formalization of the problem of parallel rearrangements calling.

The tool takes synteny blocks and a phylogenetic tree as input and outputs rearrangement events. 
The tool tests each rearrangement for consistency with a tree, and sorts the events by their parallelism score and provides diagrams of the neighbors for each block of interest, allowing the detection of horizontally transferred blocks or their extra copies and the inversions in which copied blocks are involved.

## Installation 

PaReBrick can be installed with `pip`:

```bash
pip install PaReBrick
```

Now you can run tool from any directory as `PaReBrick` (or `parebrick`).


## Script parameters
Main script of project including all modules together can be run from anywhere as console tool.

### Required input

**Important for input:** Identifiers in tree and on in blocks must be equal;

#### `--tree/t`
Tree in newick format, must be parsable by ete3 library.
You can read more about [formats supported by ete3](http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#reading-and-writing-newick-trees)

#### `--blocks_folder/-b`
Path to folder with blocks resulted as output of original Sibelia or maf2synteny tool.

[BLOCKS-OBTAIN.md](BLOCKS-OBTAIN.md) — instruction about how to obtain synteny blocks using SibeliaZ.

### Optional input

#### `--labels/-l`
Path to csv file with tree labels for showing on phylogenetic tree. 
Must contain two columns: `strain` and `label`.

#### `--output/-o`
Path to output folder.
Default is `./parebrick_output`

### Output
Output consist of tree folders.
1. `preprocessed_data` — 
here you can find both all and common synteny blocks in both `infercars`, `GRIMM` and `csv` formats.
As well as same `genomes_lengths.csv` file with lengths of provided genomes.
2. `balanced_rearrangements_output` — this folder contains `stats.csv` file with all non-convex characters statistics of balanced rearrangements. 
And folders `characters`, `tree_colorings` with character representation in rendered `.pdf` tree and `.csv` formats.
3. `unbalanced_rearrangements_output` — this folder contains `stats.csv` file with all non-convex characters statistics of unbalanced rearrangements. 
And folders `characters`, `tree_colorings` with character in rendered `.pdf` tree and `.csv` formats in subfolders according to their clustering representation.


## Example run and data
Example data is available in `example-data` folder.

### How to run example:
1. Clone all repository with data:
```bash
git clone https://github.com/ctlab/parallel-rearrangements
```

2. Change directory to data folder:
```bash
cd parallel-rearrangements/example-data/streptococcus_pyogenes/input
```

3. (a) Running tool on example input:
```
PaReBrick -t tree.nwk -b maf2synteny-output -l labels.csv
```

3. (b) Or with minimal required arguments (no labels):
```
PaReBrick -t tree.nwk -b maf2synteny-output
```
