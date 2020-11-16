## How to obtain synteny blocks using SibeliaZ

Overall idea is to use `SibeliaZ` for making alignment and then construct longer synteny blocks with `maf2synteny` module as described [here](https://github.com/medvedevgroup/SibeliaZ#building-synteny-blocks).

## 0) Data Preparation
Let's suggest you already have your genomes in `fasta` format in `merged.fna` file.

## 1) SibeliaZ Installation

We recommend installing [SibeliaZ](https://github.com/medvedevgroup/SibeliaZ) with a conda:
```bash
conda install sibeliaz
```

## 2) Running SibeliaZ Alignment
We recommend to use some parameters: 
* `-k 15` k-mer size for bacterial size genomes (recommended in SibeliaZ documentation for bacterial genomes);
* `-n` to skip the alignment with nucleotides and only output coordinates of the alignment saving time and memory;
* `-t 4` optional, for desired threads number;
* **Do not change `-m` parameter**. 
It will slow down the computation significantly and will not give you blocks. 
Blocks need to be obtained with `maf2synteny` tool.

Final command:
```bash
sibeliaz -k 15 -n -o sibeliaz_out merged.fna
```

Now you have `blocks_coords.gff` in output folder!

## 3) Constructing Synteny Blocks 

### 3.1) Preparation
We recommend using [`fine` parameters](https://github.com/bioinf/Sibelia/blob/master/SIBELIA.md#custom-parameters-set) for merging alignments into blocks.

For this, you need to create file `fine.txt` with this content:
```
30 150
100 500
500 1500
```

Also, create file `fine_500.txt` if you want to get shorter blocks (shorter then 1500):
```
30 150
100 500
````

More about these parameters and what they mean [here](https://github.com/bioinf/Sibelia/blob/master/SIBELIA.md#custom-parameters-set.

### 3.2) Running
Then, you can run `maf2syneny` SibeliaZ module with desired minimal block size `-b`. 

For getting blocks with minimal size 5000:
```
maf2synteny -s fine.txt -b 5000 blocks_coords.gff
```

For getting blocks with minimal size 1000:
```
maf2synteny -s fine_500.txt -b 1000 blocks_coords.gff
```

#### Now you have synteny blocks and can use them for PaReBrick or any of your purposes!