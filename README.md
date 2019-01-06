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

We use a specific input to illustrate how the constrained hierarchical clustering works. The raw data is derived from Hi-C experiment of RAO [ref] in GM12878 cell type using a specific 4bp-cutter restriction enzyme, MboI (SRA: SRR1658602). We selected some loci which are spread over 74Mb, 10Mb and 6Mb from the chromosome 18 at 40kb of resolution.

The data is stored inside the Master folder

```
- HTADs/Example_data/chromosome18_74Mb.tsv
- HTADs/Example_data/chromosome18_10Mb.tsv
- HTADs/Example_data/chromosome18_6Mb.tsv

```
These chromosomal regions can be represented as a symmetric contact matrixes (M), where each color point represents the average number of interactions between the bins (M[i,j]) on the chromosome. In these way, is easy to see how,as we zooming in diferents regions inside the chromosme, we can distiguish more multi-levels of TADs

![alt text](https://github.com/paulasoler/HTADs/blob/master/zoom_pictures_test-1.png)

To obtain this 40kb-binned raw interaction matrices, we processed Hi-C data using a complete Python library, called TADbit, that deals with all the necessary steps to analyze, model and explore 3C-based data. https://github.com/3DGenomes/TADbit


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
### 3.2) Pearson correlation matrix

Correlation coefficients were calculated between each bin of the matrices to construct a Pearson correlation matrix.

### 3.3) Iterative analysis
The iteraitve analysis is composed by three main steps:
- Principal component analysis
- Applying probabilistic broken-stick model
- Assesment using Calinski-Harabasz
 
A parallel analysis of the first 200 principal components (PCs) was applied. In each case, the number of clusters was extracted as the maximum between the values that are above the one pre-estimated by the probabilistic broken-stick model. For each number of clusters, the Calinski-Harabasz index was calculated. The highest value was selected, as an optimal value for dividing the data by choosing the minimum number of components that increase the dispersion among the groups generated.

#### 3.3.1) Multiples iterative analysis version

-  Fast version
NO multicore, random search para encontar la combinacion optima numero de clust i numero de PCA with arbitray precision.

- Optimazed version
Todas las combinaciones n PCA y todos los posibles numero de cluster (the best) multicore
Busqueda exhaustiva.

### 3.4) Generation of RData output

## Built With

* [R](https://www.r-project.org/about.html) - The web framework used

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Paula Soler Vila** - *Initial work* - (https://github.com/paulasoler/)
* **Pol Cuscó Pons** - *Initial work* - (https://github.com/nanakiksc/)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## References

* RAO
