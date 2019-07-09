# TADpole

TADpole, a computational tool designed to identify and analyze the entire hierarchy of topologically associated domains (TADs) in intra-chromosomal interaction matrices.

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
install.packages(c('bigmemory', 'dendextend', 'doParallel', 'foreach', 'fpc', 'Matrix', 'rioja'))
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

In this repository, we provide a publicly available HiC data set (SRA: [SRR1658602](https://www.ebi.ac.uk/ena/data/view/SRR1658602)).

In the `inst/extdata/` directory, there are 3 regions of chromosome 18 binned at 40kb, one corresponding to the full chromosome, and the others representing regions of 10 and 6 Mb:

```
- inst/extdata/chromosome18_74Mb.tsv
- inst/extdata/chromosome18_10Mb.tsv
- inst/extdata/chromosome18_6Mb.tsv
```

![Zoom](https://github.com/paulasoler/TADpole/blob/master/misc/zoom_pictures.png)

To obtain these interaction matrices, we processed the HiC data using the [TADbit](https://github.com/3DGenomes/TADbit) Python library, that deals with all the necessary steps to analyze and normalize 3C-based data.

In this tutorial, we are going to use **chromosome18_10Mb.tsv**.

### 2.1) Input data
To run the main function `TADpole`, you need to provide an intrachromosomal interaction matrix, representing an entire chromosome or a contiguous chromosome region. Input data are provided in a tab-delimited matrix format containing the interaction values in each cell. These interaction values can be the raw or normalized interaction counts. We highly recommend [ONED](https://github.com/qenvio/dryhic) normalization, as it effectively corrects for known experimental biases.


### 2.2) Running the algorithm
The basic usage is the following:

```
library(TADpole)
chromosome18_10Mb <- system.file("extdata", "chromosome18_10Mb.tsv", package = "TADpole")

tadpole <- TADpole(chromosome18_10Mb)
```

#### 2.2.1) Parameters
- **input_data**: `path` to the input file. Must be in a tab-delimited matrix format.
- **max_pcs**: the maximum number of principal components to retain for the analysis. Default value of 200 is recommended.
- **min_clusters**: minimum number of clusters into which partition the chromosome.
- **bad_frac**: fraction of the matrix to flag as bad rows/columns.
- **hist_bad_columns**: plot the distribution of row/column coverage to help in selecting a useful value for `bad_frac`. Mostly for debugging purposes.
- **centromere_search**: split the matrix by the centromere into two smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) datasets.

## 3) Output
The function `TADpole` returns a `tadpole` object containing the following items:

- ***n_pcs***: Optimal number of principal components.
- ***optimal_n_clusters***: Optimal number of clusters.
- ***dendro***: Hierarchical tree-like structure with the TAD divisions.
- ***clusters***: A list containing the TAD information of all the clusters _(x)_ defined by the broken stick model.
  + ***clusters$`x`$coord***: Start and end coordinades of the TADs.
  + ***clusters$`x`$CH-index***: Calinski-Harabasz index of this segmentation.

```
head(tadpole)

$n_pcs
[1] 35

$optimal_n_clusters
[1] 16

$dendro

Call:
rioja::chclust(d = dist(pcs))

Cluster method   : coniss
Distance         : euclidean
Number of objects: 248


$clusters
$clusters$`2`
$clusters$`2`$`CH-index`
[1] 62,12965

$clusters$`2`$coord
  start end
1     1 110
2   111 248
```

### 3.1) Plotting the results

#### 3.1.1) Dendrogram plot 
Dendrogram with all the hierarchical levels validated by the Broken-Stick model. Optimal clusters are highlighted with red rectangles. 

```
plot_dendro(tadpole)
```

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/dendogram-1_2.png" width="60%">
</p>

The optimal segmentation can be overlayed on a symmetric HiC matrix to visualize the called TADs

```
plot_borders(tadpole, chromosome18_10Mb)
```

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/TAD_partition.png" width="60%" align="center">
</p>

# DiffT Score
Difference score between topological partitions.

### 1) Input data
In the `data/` directory, there are 2 partitions from chromosome 1 obtained in two diferent conditions. Each of them
is a BED-like `data.frame`.

```
- data/control.Rdata
- data/case.Rdata
```

### 2) Computing the DiffT score

```
difft_control_case <- diffT(control, case)
```

#### 2.1) Parameters
- **bed_x**, **bed_y**: two `data.frame`s with a BED-like format with 3 columns: chromosome, start and end coordinates of each TAD, in bins.

### 3) Output
The function `diffT` returns a `numeric` vector representing the cumulative the DiffT score along the bins.
The highest local differences between the two matrices can be identified by the sharpest changes in the slope of the function.

<p align="center">
<img src="https://github.com/paulasoler/TADpole/blob/master/misc/DiffT_score.png" width="60%" align="center">
</p>

## Authors

- **Paula Soler Vila** - (https://github.com/paulasoler/)
- **Pol Cuscó Pons** - (https://github.com/nanakiksc/)

## References

1. Rao SSP, Huntley MH, Durand NC, Stamenova EK, Bochkov ID, Robinson JT, Sanborn AL, Machol I, Omer AD, Lander ES, Aiden EL. A Three-dimensional Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping. Cell. 2014;159:1665–1680.
2. Serra, F., Baù, D., Goodstadt, M., Castillo, D. Filion, G., & Marti-Renom, M.A. (2017). Automatic analysis and 3D-modelling of Hi-C data using TADbit reveals structural features of the fly chromatin colors. PLOS Comp Bio 13(7) e1005665. doi:10.1371/journal.pcbi.1005665
3. Enrique Vidal, François le Dily, Javier Quilez, Ralph Stadhouders, Yasmina Cuartero, Thomas Graf, Marc A Marti-Renom, Miguel Beato, Guillaume J Filion, OneD: increasing reproducibility of Hi-C samples with abnormal karyotypes, Nucleic Acids Research, Volume 46, Issue 8, 4 May 2018, Page e49, https://doi.org/10.1093/nar/gky064
