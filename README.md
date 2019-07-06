# TADpole

TADpole, a computational tool designed to identify and analyze the entire hierarchy of topologically associated somains (TADs) in intra-chromosomal interaction matrices.

## 1) Installation

<!--
### 1.1) Using the _devtools_ package

This is the recommended way of installing TADpole.

- First, install the devtools package from CRAN, if it is not already installed

```
install.packages("devtools")
```

- Then, Install the HTADs package from GitHub

```
devtools::install_github("paulasoler/TADpole")
```

### 1.2) Manual installation from source
-->

- First, install the required dependencies from within R

```
install.packages('bigmemory', 'Matrix','doParallel', 'dendextend', 'parallel','foreach','fpc', 'rioja')
```

- Then, get the latest version of the source code from Github

by using _wget_:

```
wget https://github.com/paulasoler/TADpole/archive/master.zip
unzip TADpole-master.zip
mv TADpole-master TADpole
```

or by cloning the repository:

```
git clone https://github.com/paulasoler/TADpole.git
```

- Finally, install the package

```
R CMD INSTALL TADpole
```

## 2) Getting started

In this tutorial, we provide a publicly available HiC data set (SRA: [SRR1658602](https://www.ebi.ac.uk/ena/data/view/SRR1658602)).

In the `data/` directory, there are 3 regions of chromosome 18 binned at 40kb, one corresponding to the full chromosome, and the others representing regions of 10Mb and 6 Mb:

```
- data/chromosome18_74Mb.tsv
- data/chromosome18_10Mb.tsv
- data/chromosome18_6Mb.tsv
```

![Zoom](https://github.com/paulasoler/TADpole/blob/master/misc/zoom_pictures.png)

To obtain this interaction matrices, we processed the HiC data using the [TADbit](https://github.com/3DGenomes/TADbit) Python library, that deals with all the necessary steps to analyze and normalize 3C-based data.

In this tutorial, we are going to use **chromosome18_10Mb.tsv**.

### 2.1) Input data
To run the main function `TADpole`, you need to provide an intrachromosomal interaction matrix, representing an entire chromosome or a contiguous chromosome region. Input data are formatted as a tab-delimited matrix format containing the interaction values in each cell. These interaction values can be the raw or normalized interaction counts. We highly recommend [ONED](https://github.com/qenvio/dryhic) normalization, as it effectively corrects for known experimental biases.


### 2.2) Running the algorithm
The basic usage is the following:
```
tadpole <- TADpole(input_data = "data/chromosome18_10Mb.tsv", max_pcs = 200, min_clusters = 2, bad_frac = 0.01, centromere_search = FALSE, hist_bad_columns = FALSE)
```

#### 2.2.1) Parameters
- **input_data**: `path` to the input file. Must be in a tab-delimited matrix format.
- **max_pcs**: Numeric. The maximum number of principal components to retain for the analysis. Default value of 200 is recommended.
- **min_clusters**: Minimum number of clusters into which partition the chromosome.
- **bad_frac**: fraction of the matrix to flag as bad rows/columns.
- **hist_bad_columns**: plot the distribution of row/column coverage to help in selecting a useful value for `bad_frac`. Mostly for debugging.
- **centromere_search**: split the matrix by the centromere into two smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) datasets.

## 3) Output
The function `TADpole` returns a `tadpole` object, which is a `list` containing the following items:

- ***n_pcs***: Optimal number of principal components.
- ***optimal_n_clusters***: Optimal number of clusters.
- ***dendro***: Hierarchical tree-like structure with the TAD divisions.
- ***clusters***: A list containing the TAD information of all the clusters _(x)_ defined by the broken stick model.
  + ***clusters$`x`$coord***: Start and end coordinades of the TADs.
  + ***clusters$`x`$CH-index***: Calinski-Harabasz index of this segmentation.

<!-- ![CHindex](https://github.com/paulasoler/HTADs/blob/master/misc/CHindex_accurate_method.png) -->

```
head(tadpole)

$n_pcs
[1] 41

$optimal_n_clusters
[1] 20

$dendro

Call:
rioja::chclust(d = dist(pcs))

Cluster method   : coniss
Distance         : euclidean
Number of objects: 251


$clusters
$clusters$`2`
$clusters$`2`$`CH-index`
[1] 62.7307

$clusters$`2`$coord
start end
1      1  26
2     27  45
3     46  72
```

### 3.1) Plotting the results

#### 3.1.1) Dendrogram plot 
Dendrogram with all the hierarchical levels validated by the Broken-Stick model.
```
plot_dendro(tadpole)
```

The optimal segmentation can be overlayed on a symmetric HiC matrix to visualize the called TADs

```
plot_borders(tadpole, input_data = "data/chromosome18_10Mb.tsv")
```

![Zoom](https://github.com/paulasoler/HTADs/blob/master/misc/dendogram-1_2.png)

## Authors

- **Paula Soler Vila** - (https://github.com/paulasoler/)
- **Pol Cuscó Pons** - (https://github.com/nanakiksc/)

## References

1. Rao SSP, Huntley MH, Durand NC, Stamenova EK, Bochkov ID, Robinson JT, Sanborn AL, Machol I, Omer AD, Lander ES, Aiden EL. A Three-dimensional Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping. Cell. 2014;159:1665–1680.
