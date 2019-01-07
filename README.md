# HTADs

The HTADs package applies a constrained hierarchical method to detect Topologically Associated Domains (TADs).

## 1) Installing of the released version of HTADs.

### 1.1) Using devtools package

- Step 1: Install the devtools package from CRAN

```
install.packages("devtools")
```
- Step 2: Install the HTADs package from GitHub

```
library(devtools)
install_github("paulasoler/HTADs")
```
- Step 3: Load the package
```
library(HTADs)
```

### 1.2) Download, unzip and install manually the package.

```
wget https://github.com/paulasoler/HTADs/archive/master.zip (you have the alternative option using clone)
unzip master.zip
mv master HTADs
sudo R CMD INSTALL HTADs
```

## 2) Dependencies

```
require(colorRamps)
require(data.table)
require(Matrix)
require(doParallel)
require(rioja)
require(fpc)
require(reshape2)
```

## 3) Getting started: Running algorithm

We use a specific input to illustrate how the constrained hierarchical clustering works. The raw data is derived from Hi-C experiment of RAO[1] in GM12878 cell type using a specific 4bp-cutter restriction enzyme, MboI (SRA: SRR1658602). We selected some loci which are spread over 74Mb, 10Mb and 6Mb from the chromosome 18 at 40kb of resolution.

The data is stored inside the Master folder

```
- HTADs/Example_data/chromosome18_74Mb.tsv
- HTADs/Example_data/chromosome18_10Mb.tsv
- HTADs/Example_data/chromosome18_6Mb.tsv
```
These chromosomal regions can be represented as a symmetric contact matrixes (M), where each color point represents the average number of interactions between the bins (M[i,j]) on the chromosome. In this way, it is easy to see (Figure 1) how, by zooming into different regions within the chromosome, we can distinguish multiple hierachical levels of TADs.

![alt text](https://github.com/paulasoler/HTADs/blob/master/zoom_pictures_test-1.png)

To obtain this 40kb-binned raw interaction matrices, we processed Hi-C data using a complete Python library, called TADbit, that deals with all the necessary steps to analyze, model and explore 3C-based data. https://github.com/3DGenomes/TADbit. 

In this tutorial, we are going to use **chromosome18_10Mb.tsv** as a template but you have the possibility to run a complete chromosome (*chromosome18_74Mb.tsv*) or a smaler data set (*chromosome18_6Mb.tsv*).


### 3.1) Format of input data:
We highly recommended normalize the data before to start the TAD caller. We give advice to use ONED normalization to correct the experimental bias that affects the signal profile defined by the contacts.
https://github.com/qenvio/dryhic

 - 3 columns-format stored in tab-separated value without header. The first and the second columns correspond to the pair of loci(i,j) which are stablish a average number of normalized interactions (the last column). 

```
28	27	1108.42577685
22	12	423.808156978
26	6	286.1442771
17	13	1740.23475628
```
### 3.2) Pearson correlation matrix and Principal Component Analysis

Correlation coefficients were calculated between each bin of the symmetric normalize matrix to construct a Pearson correlation matrix.This matrix is the input to compute the first 200 principal components*(PCs)*.

```
  # Load and clean data.
  mat <- load_mat(file_name)

  # Sparse matrix and correlation.
  sparse_matrix <- Matrix::Matrix(mat, sparse = TRUE)
  correlation_matrix <- sparse_cor(sparse_matrix)$cor
  correlation_matrix[is.na(correlation_matrix)] <- 0

  # PCA (compute first `n.pcs` components).
  number_pca <- min(max_pcs, nrow(mat))
  pca <- prcomp(correlation_matrix, rank. = number_pca)
```  

### 3.3) Find optimal clustering parameters based on Calinhara score.

#### 3.3.1) Multiples iterative analysis version

> ***Fast version***
- Based on a random permutation search that allows to find the optimal combination between the number of cluster and the number of PCs with arbitrary precision.
- Not multicore implemented.

> ***Optimazed version***
- Based of the exhaustive search method to get the best combination between all the possible PCs and number of clusters combinations.
- Multicore implemented.

![alt text](https://github.com/paulasoler/HTADs/blob/master/Optimazed_version.png)

### 3.4) Generation of RData output

The Rdata created is a R object with contains the following information:
   1) ***$n_pcs*** = Optimal number of principal components
   2) ***$optimal_n_clusters*** = Optimal number of clusters
   3) ***$dendro*** = Hierarchical tree-like strcutures with the TADs division
   4) ***$clusters$`x`$`CH-index`*** = Calinski-Harabasz index from 1 to optimal number of clusters *(x)*.
   5) ***$clusters$`x`$`coord`*** = Start and end TADs coordinades from 1 to optimal number of clusters *(x)*.


```
load("chromosome18_10Mb.Rdata")
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

## Built With

* [R](https://www.r-project.org/about.html) - The web framework used

## Authors

* **Paula Soler Vila** - *Initial work* - (https://github.com/paulasoler/)
* **Pol Cuscó Pons** - *Initial work* - (https://github.com/nanakiksc/)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## References

1. Rao SSP, Huntley MH, Durand NC, Stamenova EK, Bochkov ID, Robinson JT, Sanborn AL, Machol I, Omer AD, Lander ES, Aiden EL. A Three-dimensional Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping. Cell. 2014;159:1665–1680.

