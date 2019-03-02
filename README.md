# HTADs

The HTADs package applies a constrained hierarchical method to detect Topologically Associated Domains (TADs).

## 1) Installation

<!--
### 1.1) Using the _devtools_ package

This is the recommended way of installing HTADs.

- First, install the devtools package from CRAN, if it is not already installed

```
install.packages("devtools")
```

- Then, Install the HTADs package from GitHub

```
devtools::install_github("paulasoler/HTADs")
```

### 1.2) Manual installation from source
-->

- First, get the latest version of the source code from Github

by using _wget_:

```
wget https://github.com/paulasoler/HTADs/archive/master.zip
unzip master.zip
mv master HTADs
```

or by cloning the repository:

```
git clone https://github.com/paulasoler/HTADs.git
```

- Finally, install the package

```
R CMD INSTALL HTADs
```

## 2) Getting started

In this tutorial, we provide a publicly available HiC data set (SRA: [SRR1658602](https://www.ebi.ac.uk/ena/data/view/SRR1658602)).

In the `data/` directory, there are 3 regions of chromosome 18 binned at 40kb, one corresponding to the full chromosome, and the others representing regions of 10Mb and 6 Mb:

```
- data/chromosome18_74Mb.Rdata
- data/chromosome18_10Mb.Rdata
- data/chromosome18_6Mb.Rdata
```

![Zoom](https://github.com/paulasoler/HTADs/blob/master/misc/zoom_pictures.png)

To obtain this interaction matrices, we processed the HiC data using the [TADbit](https://github.com/3DGenomes/TADbit) Python library, that deals with all the necessary steps to analyse and normalize 3C-based datasets. 

In this tutorial, we are going to use **chromosome18_10Mb.tsv**.

### 2.1) Input data
To run the main funcion `call_HTADs`, you need to provide a `data.frame` with 3 columns. The first two columns correspond to the pair of bins _(i, j)_ and the third column is their 3C-based interaction score. This score can be the raw or the normalised interaction count. We highly recommend [ONED](https://github.com/qenvio/dryhic) normalization, as it effectively corrects for known experimental biases. 

```
28 27 1108.4257768
22 12 423.8081569
26 6  286.1442771
17 13 1740.2347562
```

### 2.2) Running the algorithm
The basic usage is the following:
```
load('data/chromosome18_10Mb.Rdata')
htads <- call_HTADs(chromosome18_10Mb, method = 'accurate')
```

#### 2.2.1) Parameters
- **input_data**: `data.frame` with 3 columns containing the 3C-based dataset in the format `(bin1, bin2, score)`.
- **cores**: Numeric. When `method` is `"accurate"`, the number of cores to use for parallel execution.
- **max_pcs**: Numeric. The maximum number of principal components to retain for the analysis. Default value of 200 is recommended.
- **method**: Character. Which version of the algorithm to use. `"accurate"` (default) or `"fast"`.

  -- ***Accurate***
    - Performs an _exhaustive_ search of all possible solutions.
    - Has support for multicore parallelization.

  -- ***Fast***
    - Performs a random search over all possible solutions.
    - Its precision can be increased at the expense of computation time by increasing the `n_samples` parameter.

- **n_samples**: Numeric. When `method` is `"fast"`, the number of samples used to approximate the optimal solution.
- **min_clusters**: Minimum number of clusters into which partition the chromosome.
- **plot**: Logical. Whether to plot the scores of every tested `n_pcs`/`n_clusters` combination.

## 3) Output
The funcion `call_HTADs` returns a `htads` object, which is a `list` contaning the following items:

- ***n_pcs***: Optimal number of principal components.
- ***optimal_n_clusters***: Optimal number of clusters.
- ***dendro***: Hierarchical tree-like structure with the TAD divisions.
- ***clusters***: A list containing the TAD information of all the clusters _(x)_ defined by the broken stick model.
  + ***clusters$`x`$coord***: Start and end coordinades of the TADs.
  + ***clusters$`x`$CH-index***: Calinski-Harabasz index of this segmentation.

<!-- ![CHindex](https://github.com/paulasoler/HTADs/blob/master/misc/CHindex_accurate_method.png) -->

```
head(htads)

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
The optimal segmentation can be overlayed on a symmetric HiC matrix to visualize the called TADs

```
plot_borders(matrix_chr18_10Mb, htads)
```

![Zoom](https://github.com/paulasoler/HTADs/blob/master/misc/dendogram-1_2.png)

## Authors

- **Paula Soler Vila** - (https://github.com/paulasoler/)
- **Pol Cuscó Pons** - (https://github.com/nanakiksc/)

## References

1. Rao SSP, Huntley MH, Durand NC, Stamenova EK, Bochkov ID, Robinson JT, Sanborn AL, Machol I, Omer AD, Lander ES, Aiden EL. A Three-dimensional Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping. Cell. 2014;159:1665–1680.
