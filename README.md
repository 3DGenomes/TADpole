# TADpole

TADpole is a computational tool designed to identify and analyze the entire hierarchy of topologically associated domains (TADs) in intra-chromosomal interaction matrices.

## 1) Installation

<!--
### 1.1) Using the _devtools_ package

This is the recommended installation procedure.

- First, install the devtools package from CRAN

```R
install.packages("devtools")
```

- Then, Install the TADpole package from GitHub

```R
devtools::install_github("paulasoler/TADpole")
```

### 1.2) Manual installation from source
-->

- First, install the required dependencies in R

```R
install.packages(c('bigmemory', 'cowplot', 'doParallel', 'foreach', 'fpc',
                   'ggdendro', 'ggplot2', 'ggpubr', 'gridExtra', 'Matrix',
                   'plyr', 'reshape2', 'rioja', 'viridis', 'zoo'))
```

- Then, get the latest version of the source code from Github

by using _wget_:

```Shell
wget https://github.com/paulasoler/TADpole/archive/master.zip
unzip master.zip
mv TADpole-master TADpole
```

or by cloning the repository:

```Shell
git clone https://github.com/paulasoler/TADpole.git
```

- Finally, install TADpole.

```Shell
R CMD INSTALL TADpole
```

Note: if you download the zip file from the GitHub website instead, it will be named `TADpole-master`, so please adapt the `unzip` command accordingly.

## 2) Getting started

In this repository, we provide a test case from a publicly available Hi-C data set (SRA: [SRR1658572](https://www.ebi.ac.uk/ena/data/view/SRR1658572)) (1).

In the `inst/extdata/` directory, we provided a 3 Mb region (chr18:9,200,000-12,200,000) of a human Hi-C dataset at 30kb resolution. 

```
inst/extdata/raw_chr18_300_500_30kb.tsv
```

To obtain this interaction matrix, we processed the Hi-C data using the [TADbit](https://github.com/3DGenomes/TADbit) (2) Python package, that deals with all the necessary processing and normalization steps.

### 2.1) Input data
To run the main function `TADpole`, you need to provide an intrachromosomal interaction matrix, representing an entire chromosome or a chromosome region. The input is a tab-separated values file containing the interaction matrix (_M_) with _N_ rows and _N_ columns, where _N_ is the number of bins in which the chromosome region is divided. Each position of the matrix (_M<sub>ij</sub>_) contains the interaction values (raw or normalized) between the corresponding pair of genomic bins _i_ and _j_. We recommend [ONED](https://github.com/qenvio/dryhic) (3) normalization, as it effectively corrects for known experimental biases.

### 2.2) Running the algorithm

Schematic overview of the TADpole algorithm (for further details, refer to [Soler-Vila _et al_.](https://www.biorxiv.org/content/10.1101/698720v1) (4)

![Zoom](https://github.com/paulasoler/TADpole/blob/master/misc/Figure1.png)

The basic usage is the following:

```R
library(TADpole)
mat_file <- system.file("extdata", "raw_chr18_300_500_30kb.tsv", package = "TADpole")

tadpole <- TADpole(mat_file, chr = "chr18", start = 9000000, end = 15000000, resol = 30000)
```

#### 2.2.1) Parameters
- **mat_file**: path to the input file. The file must be in a tab-delimited matrix format.
- **chr**: chromosome name.
- **start**: initial position of the chromosomal region or the chromosome, in base pairs.
- **end**: final position of the chromosomal region or the chromosome, in base pairs.
- **resol**: binning-size of the Hi-C experiment, in base pairs.
- **max_pcs**: the maximum number of principal components to retain for the analysis. The default (recommended) value is 200.
- **min_clusters**: minimum number of TADs to find.
- **bad_frac**: fraction of the matrix to flag as bad columns.
<!--
- **hist_bad_columns**: `logical`. Plot the distribution of column coverage to help in selecting a useful value for `bad_frac`. Mostly for debugging purposes.
-->
- **centromere_search**: `logical`. Split the matrix by the centromere into two sub-matrices representing the chromosome arms. Useful when working with big matrices (>15000 bins).

## 3) Output
The function `TADpole` returns a `tadpole` object containing the following descriptors:

- **n_pcs**: optimal number of principal components.
- **optimal_n_clusters**: optimal number of chromatin partitions (that is the index of the optimal level + 1).
- **dendro**: hierarchical tree-like structure cut at the maximum significant number of levels identified by the broken-stick model (max(ND)).
- **clusters**: a `list` containing the TADs for each hierarchical level _(x)_ defined by the broken stick model.
  + **clusters$`x`**: start and end coordinades of all TADs.
- **score**: CH index associated to each dendrogram.
- **merging_arms**: if `centromere_search` is `TRUE`, contains the start and end coordinates of the TADs of the full chromosome.

```R
head(tadpole)

$n_pcs
[1] 20

$optimal_n_clusters
[1] 12

$dendro

Call:
rioja::chclust(d = dist(pcs))

Cluster method   : coniss
Distance         : euclidean
Number of objects: 198

$clusters
$clusters$`2`
  start end
1     1 110
2   111 200
...

$scores
     1        2        3        4        5        6        7        8        9
1   NA 47,90916 42,22857 39,40353 43,61547 41,24569  0,00000  0,00000  0,00000
2   NA 44,47879 43,28183 45,06219 44,02830 45,38542 49,09032  0,00000  0,00000
...

```

### 3.1) Plotting the results

#### 3.1.1) Raw Hi-C map and histogram of the interaction values.
Automatically, TADpole generates a map of the intra-chromosomal interaction matrix under study, together with a histogram showing the distribution of interaction values. In the latter, a dashed line that indicates the number of columns (and corresponding rows) excluded from the analysis for having a low number of interactions (the so-called bad columns). Specifically, the columns (and rows) that contain an empty cell at the main diagonal, and those whose cumulative interactions are below the first (by default) percentile, are excluded from the analysis.

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/Figure2.png" width="70%">
</p>

#### 3.1.2) Hierarchical plot
**Left**, the complete dendrogram obtained from the Hi-C matrix cut at a maximum significant number of levels (_max(ND)_) reported by the broken-stick model (including the partitions in 2 up to 16 TADs). Among these levels, the highest-scoring one is selected according to the CH index analysis. **Right**, Hi-C contact map showing the complete hierarchy of the significant levels selected by the broken stick model (black lines) along with the optimal one with 12 TADs, identified by the highest CH index (blue line).

```R
plot_hierarchy(mat_file, tadpole, chr = "chr18", start = 9000000, end = 15000000, resol = 30000)
```

##### 3.1.2.1) Parameters
- **mat_file**: path to the input file. The file must be in a tab-delimited matrix format.
- **tadpole**: `tadpole` object
- **chr**: chromosome name.
- **start**: initial position of the chromosomal region or the chromosome, in base pairs.
- **end**: final position of the chromosomal region or the chromosome, in base pairs.
- **resol**: binning-size of the Hi-C experiment, in base pairs.
- **centromere_search**: `logical`. Split the matrix by the centromere into two sub-matrices representing the chromosome arms. Useful when working with big matrices (>15000 bins).

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/Figure3_1.png" width="70%">
</p>

#### 3.1.3) Matrix of Calinski-Harabasz indices

```R
CH_map(tadpole)
```

##### 3.1.3.1) Parameters
- **tadpole**: `tadpole` object.

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/Figure4.png" width="50%" align="center">
</p>

# DiffT score
To compare pairs of topological partitions, _P_ and _Q_, identified by TADpole at the same level of the hierarchy, we defined a Difference Topology score (DiffT). Specifically, the partitioned matrices are transformed into binary forms _p_ for _P_, and _q_ for _Q_, in which each entry _p<sub>ij</sub>_ (_q<sub>ij</sub>_) is equal to 1 if the bins _i_ and _j_ are in the same TAD and 0 otherwise. Then, the DiffT is computed as the normalized (from 0 to 1) difference between the binarized matrices as a function of the bin index _b_ as:

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/DiffT_formula.png" width="30%" align="center">
</p>

where _N_ is the total number of bins.
<br>

### 1) Input data
Here, the DiffT score analysis is used to compare the chromatin partitions at the same hierarchical level determined in two different experiments: control and case.

In the `inst/extdata/` directory, there are 2 files in a BED-like format.

```
inst/extdata/control.bed
inst/extdata/case.bed
```

### 2) Computing the DiffT score

```R
control <- read.table(system.file("extdata", "control.bed", package = "TADpole"))
case <- read.table(system.file("extdata", "case.bed", package = "TADpole"))

difft_control_case <- diffT(control, case)
```

#### 2.1) Parameters
- **bed_x**, **bed_y**: two `data.frame`s with a BED-like format with 3 columns: chromosome, start and end coordinates of each TAD, in bins.

### 3) Output
The function `diffT` returns a `numeric` vector representing the cumulative DiffT score profiles as a function of the matrix bins.
The highest local differences between the two matrices can be identified by the sharpest changes in the slope of the function.

```R
plot(difft_control_case, type = "l")
``````

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/DiffT_score.png" width="60%" align="center">
</p>

## Authors

- **Paula Soler Vila** - ([@paulasoler](https://github.com/paulasoler))
- **Pol Cuscó Pons** - ([@nanakiksc](https://github.com/nanakiksc))
- **Marco Di Stefano** - ([@MarcoDiS](https://github.com/MarcoDiS))

## References

1. RAO, Suhas SP, et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell, 2014, 159.7: 1665-1680.
2. SERRA, François, et al. Automatic analysis and 3D-modelling of Hi-C data using TADbit reveals structural features of the fly chromatin colors. PLoS computational biology, 2017, 13.7: e1005665.
3. VIDAL, Enrique, et al. OneD: increasing reproducibility of Hi-C samples with abnormal karyotypes. Nucleic acids research, 2018, 46.8: e49-e49.
4. SOLER-VILA, Paula, et al. Hierarchical chromatin organization detected by TADpole. biorXiv, doi: https://doi.org/10.1101/698720
