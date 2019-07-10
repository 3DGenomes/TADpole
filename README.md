# TADpole

TADpole is a computational tool designed to identify and analyze the entire hierarchy of topologically associated domains (TADs) in intra-chromosomal interaction matrices.

## Authors
TADpole is currently developed at the [MarciusLab](http://www.marciuslab.org) with the contributions of [Paula Soler](https://github.com/paulasoler/), [Pol Cuscó](https://github.com/nanakiksc/), [Marco Di Stefano](https://github.com/MarcoDis/), [Irene Farabella](https://github.com/iosonoirene/), and many other members of our Lab.

## 1) Installation

<!--
### 1.1) Using the _devtools_ package

This is the recommended installation procedure.

- First, install the devtools package from CRAN

```
install.packages("devtools")
```

- Then, Install the TADpole package from GitHub

```
devtools::install_github("paulasoler/TADpole")
```

### 1.2) Manual installation from source
-->

- First, install the required dependencies in R

```
install.packages(c('bigmemory', 'dendextend', 'doParallel', 'foreach', 'fpc', 'Matrix', 'rioja'))
```

- Then, get the latest version of the source code from Github

by using _wget_:

```
wget https://github.com/3DGenomes/TADpole/archive/master.zip
unzip TADpole-master.zip
mv TADpole-master TADpole
```

or by cloning the repository:

```
git clone https://github.com/3DGenomes/TADpole.git
```

- Finally, install the TADpole package. The package 'TADpole' requires R >= 3.5.2


```
R CMD INSTALL TADpole
```

## 2) Getting started

In this repository, we provide a test case from a publicly available Hi-C data set (SRA: [SRR1658602](https://www.ebi.ac.uk/ena/data/view/SRR1658602)) (1).

In the `inst/extdata/` directory, there are 3 regions of chromosome 18 binned at 40kb, one corresponding to the full chromosome, and the others representing regions of 10 and 6 Mb:

```
- inst/extdata/chromosome18_74Mb.tsv
- inst/extdata/chromosome18_10Mb.tsv
- inst/extdata/chromosome18_6Mb.tsv
```

![Zoom](https://github.com/3DGenomes/TADpole/tree/master/misc/zoom_pictures.png)


To obtain these interaction matrices, we processed the Hi-C data using the [TADbit](https://github.com/3DGenomes/TADbit) (2) Python library, that deals with all the necessary steps to analyze and normalize Hi-C data.

In this tutorial, we are going to use **chromosome18_10Mb.tsv**.

### 2.1) Input data
To run the main function `TADpole`, you need to provide an intrachromosomal interaction matrix, representing an entire chromosome or a contiguous chromosome region. Input data are provided in a tab-delimited matrix format containing the interaction values in each cell. These interaction values can be the raw or normalized interaction counts. We recommend [ONED](https://github.com/qenvio/dryhic) (3) normalization, as it effectively corrects for known experimental biases.


### 2.2) Running the algorithm
The basic usage is the following:

```
library(TADpole)
chromosome18_10Mb <- system.file("extdata", "chromosome18_10Mb.tsv", package = "TADpole")

tadpole <- TADpole(chromosome18_10Mb)
```

#### 2.2.1) Parameters
- **input_data**: `path` to the input file. Must be in a tab-delimited matrix format.
- **max_pcs**: `numeric` the maximum number of principal components to retain for the analysis. Default value of 200 is recommended.
- **min_clusters**: `numeric` minimum number of clusters to partition the chromatin region.
- **bad_frac**: `numeric` fraction of the matrix to flag as bad columns.
- **hist_bad_columns**: `logical` plot the distribution of column coverage to help in selecting a useful value for `bad_frac`. Mostly for debugging purposes.
- **centromere_search**: `logical` split the matrix by the centromere into two smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) matrices.

## 3) Output
The function `TADpole` returns a `tadpole` object containing the following items:

- ***n_pcs***: optimal number of principal components.
- ***optimal_n_clusters***: optimal number of clusters.
- ***dendro***: hierarchical tree-like structure with the TAD divisions.
- ***clusters***: a list containing the TAD information of all the clusters _(x)_ defined by the broken stick model.
  + ***clusters$`x`***: start and end coordinades of the TADs.
- ***merging_arms***: if `centromere_search` is `TRUE`, contains the start and end coordinates of the TADs of the full chromosome.

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
  start end
1     1 110
2   111 248
```

### 3.1) Plotting the results

#### 3.1.1) Dendrogram plot
Dendrogram with all the hierarchical levels validated by the Broken-Stick model. The optimal clusters are highlighted with red rectangles.

```
plot_dendro(tadpole)
```
##### 3.1.1.1) Parameters
- **tadpole**: `tadpole` object
- **centromere_search**: `logical` split the matrix by the centromere into two smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) matrices.

<p align="center">
<img src="https://github.com/3DGenomes/TADpole/tree/master/misc/dendogram-1_2.png" width="30%">
</p>



#### 3.1.2) Optimal segmentation overlayed on the Hi-C matrix

```
plot_borders(tadpole, chromosome18_10Mb)
```
##### 3.1.2.1) Parameters
- **tadpole**: `tadpole` object
- **input_data**: `path` to the input file. Must be in a tab-delimited matrix format.
- **centromere_search**: `logical` split the matrix by the centromere into two smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) matrices.


<p align="center">
<img src="https://github.com/3DGenomes/TADpole/tree/master/misc/TAD_partition.png" width="60%" align="center">
</p>

# DiffT Score
Difference score between topological partitions.

### 1) Input data
In the `data/` directory, there are 2 partitions from mouse chromosome 1 obtained in two different conditions (4). Each of them is a BED-like `data.frame`.

```
- data/control.bed
- data/case.bed
```

### 2) Computing the DiffT score

```
control <- read.table(system.file("extdata", "control.bed", package = "TADpole"))
case <- read.table(system.file("extdata", "case.bed", package = "TADpole"))

difft_control_case <- diffT(control, case)
```

#### 2.1) Parameters
- **bed_x**, **bed_y**: two `data.frame`s with a BED-like format with 3 columns: chromosome, start and end coordinates of each TAD, in bins.

### 3) Output
The function `diffT` returns a `numeric` vector representing the cumulative the DiffT score along the bins.
The highest local differences between the two matrices can be identified by the sharpest changes in the slope of the function.
```
plot(difft_control_case, type="l")
```
<p align="center">
<img src="https://github.com/3DGenomes/TADpole/tree/master/misc/DiffT_score.png" width="60%" align="center">
</p>

## Citation

Please, cite this article if you use TADpole.

Soler-Vila, P., Cuscó, P., Farabella, I., Di Stefano, M., & Marti-Renom, M.A. (2019).
**Hierarchical levels of chromatin organization detected by TADpole**
(in preparation)

## References

1. RAO, Suhas SP, et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. Cell, 2014, 159.7: 1665-1680.
2. SERRA, François, et al. Automatic analysis and 3D-modelling of Hi-C data using TADbit reveals structural features of the fly chromatin colors. PLoS computational biology, 2017, 13.7: e1005665.
3. VIDAL, Enrique, et al. OneD: increasing reproducibility of Hi-C samples with abnormal karyotypes. Nucleic acids research, 2018, 46.8: e49-e49.
4. KRAFT, Katerina, et al. Serial genomic inversions induce tissue-specific architectural stripes, gene misexpression and congenital malformations. Nature cell biology, 2019, 21.3: 305.
