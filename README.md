# Parallel rearrangements

## Project goals
Bacterial genomes are remarkably plastic on the evolutionary time scale, and genomic rearrangements such as inversions, deletions, insertions, and duplications independently occur in the genomes of different strains, which may serve a mechanism of the adaptation to changing environmental conditions. 
Identification of these events requires laborious manual inspection and verification of the phyletic pattern consistency. 
Thus, the main goal of this project is identification and description of parallel rearrangements occurring in bacterial genomes.


## Methods
This method takes synteny blocks and phylogenetic tree for some stamps as input and 
Method consist of two main parts:
1. Constructing characters and them states (~colors of leaves):
    * Balanced rearrangements (focused in inversions) — search pattern in multiple breakpoint graph;
    * Unbalanced rearrangements — considered as copy number variation.
2. Checking convexity of constructed characters on a phylogenetic tree. 
Also counting "measure" of non-convexity for all non-convex characters and sort by its value.

For output size reduction and detection of phyletic patterns of blocks in unbalanced rearrangements case, character clustering is performed.

## System requirements
For running this project code you need to have `Python 3` with those libraries installed:
`bg`, `sklearn`, `ete3`, `numpy`, `pandas`.

## Script parameters
Main script of project including all modules together can be run with python and console arguments.
### Input
Here is input parameters of that script:
```bash
PYTHONPATH=. python3 src/main.py --help
usage: main.py [-h] --tree TREE --blocks_folder BLOCKS_FOLDER
               [--labels LABELS] --output OUTPUT

Based on synteny blocks and phylogenetic tree this tool calls parallel
rearrangements.

optional arguments:
  -h, --help            show this help message and exit
  --tree TREE, -t TREE  Tree in newick format, must be parsable by ete3
                        library.You can read more about formats supported by
                        ete3 at http://etetoolkit.org/docs/latest/tutorial/tut
                        orial_trees.html#reading-and-writing-newick-trees
  --blocks_folder BLOCKS_FOLDER, -b BLOCKS_FOLDER
                        Path to folder with blocks resulted as output of
                        original Sibelia or maf2synteny tool.
  --labels LABELS, -l LABELS
                        Path to csv file with tree labels, must contain two
                        columns: `strain` and `label`.
  --output OUTPUT, -o OUTPUT
                        Path to output folder.
```

### Output
Output consist of tree folders.
1. `preprocessed_data` — 
here you can find `blocks_coords.infercars` for better blocks representation and region annotation.
And `genomes_permutations_unique.txt` file used for calling balanced characters.
2. `balanced_rearrangements_output` — this folder contains `stats.csv` file with all non-convex characters statistics of balanced rearrangements. 
And folders `characters`, `trees` with character representation in rendered `.pdf` tree and `.csv` formats.
3. `unbalanced_rearrangements_output` — this folder contains `stats.csv` file with all non-convex characters statistics of unbalanced rearrangements. 
And folders `characters`, `trees` with character in rendered `.pdf` tree and `.csv` formats in subfolders according to their clustering representation.


## Example
Code must be with `PYTHONPATH` in root of the project. 
Running on example input:
```
PYTHONPATH=. python3 src/main.py -t example-data/e_coli/input/tree.nwk -b example-data/e_coli/input/maf2synteny-output -o example-data/e_coli/output -l example-data/e_coli/input/labels.csv
```

#### Example output
```
<Running module: Preprocess Data>
-------------------------------------------------------------------------------- 
Converting block coords to infercars format
-------------------------------------------------------------------------------- 
<Module finished: Preprocess Data>, elapsed 8.155789852142334 seconds


<Running module: Parsers and check strains>
-------------------------------------------------------------------------------- 
Strains in blocks file count: 414
Strains in tree leafs count: 414
Intersected strains count: 414
-------------------------------------------------------------------------------- 
<Module finished: Parsers and check strains>, elapsed 0.5608758926391602 seconds


<Running module: Balanced rearrangements characters>
-------------------------------------------------------------------------------- 
Breakpoint graph parsed
Getting characters from breakpoint graph component, size=66
-------------------------------------------------------------------------------- 
<Module finished: Balanced rearrangements characters>, elapsed 27.029563188552856 seconds


<Running module: Balanced rearrangements stats>
-------------------------------------------------------------------------------- 
Counting distances between unique one-copy blocks, may take a while
Got characters after breakpoint graph consideration: 33
Left non-convex characters after filtering: 26
-------------------------------------------------------------------------------- 
<Module finished: Balanced rearrangements stats>, elapsed 7.653030157089233 seconds


<Running module: Balanced rearrangements output>
-------------------------------------------------------------------------------- 
-------------------------------------------------------------------------------- 
<Module finished: Balanced rearrangements output>, elapsed 27.654157876968384 seconds


<Running module: Unbalanced rearrangements characters>
-------------------------------------------------------------------------------- 
-------------------------------------------------------------------------------- 
<Module finished: Unbalanced rearrangements characters>, elapsed 0.05320596694946289 seconds


<Running module: Unbalanced rearrangements stats and clustering>
-------------------------------------------------------------------------------- 
Got characters after copy number variation consideration: 409
Left non-convex characters after filtering: 261
Counting distances between non-convex character blocks, may take a while
Jaccard index matrix constructed
Proximity matrix constructed
Clustring is done
Clusters: 86
-------------------------------------------------------------------------------- 
<Module finished: Unbalanced rearrangements stats and clustering>, elapsed 16.02768898010254 seconds


<Running module: Unbalanced rearrangements output>
-------------------------------------------------------------------------------- 
-------------------------------------------------------------------------------- 
<Module finished: Unbalanced rearrangements output>, elapsed 214.91788911819458 seconds


Total elapsed time: 302.05323600769043
```
