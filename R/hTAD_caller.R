#' Population adjuster
#'
#' This function appends state population data
#' @param file_name The path to the file to read.
#' @param cores The number of cores to use for parallel execution.
#' @param n_pcs The maximum number of principal components to retain for the analysis.
#' @param random Logical. Whether to use the fast version of the algorithm.
#' @param n_samples When random is TRUE, the number of samples used to approximate the optimal solution.ç
#' @param plot Logical. Whether to plot the
#' @keywords per capita
#' @import dplyr
#' @export
#' @examples
#' call_hTADs("file_name.abc")

# TODO: either create one file per function (nice). Or one docstring just before each function.
# I understand only exported functions in NAMESPACE are to have docs (maybe only call_hTADs).

# TODO: refer to external functions as package::function(args)
# The :: operator adds a small overhead.

# Load packages.
require(data.table)
require(dryhic)
require(Matrix)
require(doParallel)
require(rioja)
require(fpc)
require(colorRamps)

sparse_cor <- function(x) {
  # Create a sparse correlation matrix.
  covmat <- (as.matrix(crossprod(x)) - nrow(x) * tcrossprod(colMeans(x))) / (nrow(x) - 1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat / tcrossprod(sdvec)
  list(cov = covmat, cor = cormat)
}

find_params_accurate <- function(pca, number_pca, cores) {
  registerDoParallel(cores = cores)
  calinhara_score <- foreach(i = 1:number_pca) %dopar% {
    pcs <- as.matrix(pca$x[, 1:i])
    row.names(pcs) <- 1:nrow(pcs)

    clust <- rioja::chclust(dist(pcs))

    # Broken stick.
    bs <- rioja::bstick(clust, nrow(pcs) - 1, plot = FALSE)
    r <- rle(bs$dispersion > bs$bstick)
    n_cluster <- r$lengths[r$values][1]
    chc <- cutree(clust, k = n_cluster)

    score <- rep(NA, n_cluster)
    for (n in 2:n_cluster) {
      chc <- cutree(clust, k = n)
      score[n] <- fpc::calinhara(pca$x, chc, cn = n)
    }

    score
  }

  scores <- matrix(NA, nrow = length(calinhara_score), ncol = max(sapply(calinhara_score, length)))
  for (pc in 1:length(calinhara_score)) scores[pc, 1:length(calinhara_score[[pc]])] <- calinhara_score[[pc]]
  rownames(scores) <- 1:number_pca
  colnames(scores) <- 1:ncol(scores) + 1

  optimal_PCs <- which.max(rowMeans(scores, na.rm = TRUE))
  optimal_n_clusters <- which.max(scores[optimal_PCs, ])
  message(paste('Optimal number of PCs:', optimal_PCs))
  message(paste('Optimal number of clusters:', optimal_n_clusters))

  list(n_PCs = optimal_PCs, n_clusters = optimal_n_clusters, scores = scores)
}

find_params_fast <- function(pca, number_pca, cores, n_samples) {
  n_pca <- sample(number_pca, n_samples, replace = TRUE)
  n_clu <- c()

  registerDoParallel(cores = cores)
  calinhara_score <- foreach(i = 1:n_samples) %dopar% {
    pcs <- as.matrix(pca$x[, 1:n_pca[i]])
    row.names(pcs) <- 1:nrow(pcs)

    clust <- rioja::chclust(dist(pcs))

    # Broken stick.
    bs <- rioja::bstick(clust, nrow(pcs) - 1, plot = FALSE)
    r <- rle(bs$dispersion > bs$bstick)
    n_cluster <- r$lengths[r$values][1]
    chc <- cutree(clust, k = n_cluster)

    n_clu <<- c(n_clu, sample(2:n_cluster, 1))
    chc <- cutree(clust, k = tail(n_clu, 1))
    fpc::calinhara(pca$x, chc, cn = tail(n_clu, 1))
  }

  scores <- matrix(NA, nrow = max(n_pca), ncol = max(n_clu))
  for (i in 1:n_samples) scores[n_pca[i], n_clu[i]] <- calinhara_score[[i]]
  rownames(scores) <- 1:max(n_pca)
  colnames(scores) <- 1:max(n_clu)

  optimal_PCs <- which.max(rowMeans(scores, na.rm = TRUE))
  optimal_n_clusters <- which.max(scores[optimal_PCs, ])
  message(paste('Optimal number of PCs:', optimal_PCs))
  message(paste('Optimal number of clusters:', optimal_n_clusters))

  list(n_PCs = optimal_PCs, n_clusters = optimal_n_clusters, scores = scores)
}

plot_scores <- function(optimal_params) {
  # Plot nPCs vs nClusters CHi.
  s <- optimal_params$scores
  image(s, col = green2red(1024), xaxt = 'n', yaxt = 'n', xlab = '# PCs', ylab = '# clusters')
  axis(1, at = seq(0, 1, length.out = nrow(s)), labels = rownames(s))
  axis(2, at = seq(0, 1, length.out = ncol(s)), labels = colnames(s))
}

call_hTADs <- function(file_name, cores = 1, n_pcs = 200, random = FALSE, n_samples = 60, plot = FALSE) {
  # Load and clean data.
  # mat <- as.matrix(fread(file_name, drop = 1:2)) # TODO: Drop first 2 columns?
  mat <- as.matrix(fread(file_name))
  mat <- mat[2:293, 2:293]
  mat <- symmetrize_matrix(mat)
  mat[is.na(mat)] <- 0 # Clean NA/NaN values.

  # Sparse matrix and correlation.
  sparse_matrix <- Matrix(mat, sparse = TRUE)
  correlation_matrix <- sparse_cor(sparse_matrix)$cor
  correlation_matrix[is.na(correlation_matrix)] <- 0 # Clean NA/NaN values.

  # PCA (compute first `n.pcs` components).
  number_pca <- min(n_pcs, nrow(mat))
  pca <- prcomp(correlation_matrix, rank. = number_pca)

  # Find optimal clustering parameters based on Calinhara score.
  if (random) optimal_params <- find_params_fast(pca, number_pca, cores, n_samples)
  else optimal_params <- find_params_accurate(pca, number_pca, cores)

  # Plot nPCs vs nClusters Calinhara score.
  if (plot) plot_scores(optimal_params)

  # Cluster the PCs subset with the best mean-CHI criterion.
  pcs <- as.matrix(pca$x[, 1:optimal_params$n_PCs])
  row.names(pcs) <- 1:nrow(pcs)
  clust <- rioja::chclust(dist(pcs))

  # Create hierarchical cluster object.
  htads <- list('n_pcs' = optimal_params$n_PCs[[1]],
                'optimal_n_clusters' = optimal_params$n_clusters[[1]],
                'dendro' = clust,
                'clusters' = list())
  for (n in which(!is.na(optimal_params$scores[optimal_params$n_PCs, ]))) {
    eb <- cumsum(table(cutree(clust, k = n)))
    htads$clusters[[as.character(n)]] <- list('CH-index' = optimal_params$scores[optimal_params$n_PCs, n],
                                              'coord' = data.frame('start' = c(1, eb[-length(eb)] + 1, use.names = FALSE),
                                                                   'end' = eb))
  }

  htads
}
