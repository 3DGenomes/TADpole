#' Convert a 3-column data frame to a symmetric HiC matrix
#'
#' @param input_data `data.frame` with 3 columns containing HiC data in the format `(bin1, bin2, score)`.
#' @export
load_mat <- function(input_data){
    colnames(input_data) <- paste0('V', 1:3)
    mat <- reshape2::acast(input_data, V1 ~ V2, value.var = 'V3')
    mat[is.na(mat)] <- 0 # Clean NA/NaN values.
    as.matrix(Matrix::forceSymmetric(mat, uplo = 'L'))
}

sparse_cor <- function(x) {
  # Create a sparse correlation matrix.
  covmat <- (crossprod(as.matrix(x)) - nrow(x) * tcrossprod(colMeans(x))) / (nrow(x) - 1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat / tcrossprod(sdvec)
  list(cov = covmat, cor = cormat)
}

find_params_accurate <- function(pca, number_pca, cores, min_clusters) {
  doParallel::registerDoParallel(cores = cores)
  calinhara_score <- foreach::foreach(i = 1:number_pca) %dopar% {
    pcs <- as.matrix(pca$x[, 1:i])
    row.names(pcs) <- 1:nrow(pcs)

    clust <- rioja::chclust(dist(pcs))

    # Broken stick.
    bs <- rioja::bstick(clust, nrow(pcs) - 1, plot = FALSE)
    r <- rle(bs$dispersion > bs$bstick)
    n_cluster <- r$lengths[r$values][1]

    score <- rep(NA, n_cluster)
    for (n in min_clusters:n_cluster) {
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

find_params_fast <- function(pca, number_pca, n_samples, min_clusters) {
  n_pca <- sample(number_pca, n_samples, replace = TRUE)
  n_clu <- rep(NA, n_samples)
  calinhara_score <- list()

  for (i in 1:n_samples) {
    pcs <- as.matrix(pca$x[, 1:n_pca[i]])
    row.names(pcs) <- 1:nrow(pcs)

    clust <- rioja::chclust(dist(pcs))

    # Broken stick.
    bs <- rioja::bstick(clust, nrow(pcs) - 1, plot = FALSE)
    r <- rle(bs$dispersion > bs$bstick)
    n_cluster <- r$lengths[r$values][1]

    n_clu[i] <- sample(min_clusters:n_cluster, 1)
    chc <- cutree(clust, k = n_clu[i])
    calinhara_score[[i]] <- fpc::calinhara(pca$x, chc, cn = n_clu[i])
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

#' Plot a heatmap of the Caliski-Harabasz index of every tested `n_pcs`/`n_clusters` combination
#'
#' @param htads `htad` object returned by `call_HTADs`.
#' @param file path to save the plot in html format. If `NULL` (default), the plot will not be saved on disk.
#' @examples
#' load('data/chromosome18_10Mb.Rdata')
#' htads <- call_HTADs(chromosome18_10Mb)
#' plot_scores(htads, file = 'scores_chromosome18_10Mb.html')
#' @export

plot_scores <- function(htads, file = NULL) {
  # Plot nPCs vs nClusters CHi.
  s <- htads$scores
  p <- plotly::layout(plotly::plot_ly(z = s, type = "heatmap"),
                      title = 'Caliski-Harabasz index',
                      xaxis = list(title = 'Number of clusters'),
                      yaxis = list(title = 'Number of PCs'))
  if (!is.null(file)) htmlwidgets::saveWidget(p, file)
  p
}

#' Plot TAD borders on a HiC matrix
#'
#' @param htads `htad` object returned by `call_HTADs`.
#' @param input_data `data.frame` with 3 columns containing HiC data in the format `(bin1, bin2, score)`.
#' @param add.dendro Logical. If `TRUE`, plot the clustering dendrogram beside the HiC matrix.
#' @examples
#' load('data/chromosome18_10Mb.Rdata')
#' htads <- call_HTADs(chromosome18_10Mb)
#' plot_borders(htads, chromosome18_10Mb)
#' @export

plot_borders <- function(htads, input_data) {
    mat <- load_mat(input_data)
    start_coord <- htads$clusters[[as.character(htads$optimal_n_clusters)]]$coord$start
    colors <- colorRampPalette(c('white', 'firebrick3'))

    lattice::levelplot(as.matrix(log(mat)), col.regions = colors, scales = list(draw = FALSE), colorkey = FALSE,
                       xlab = NULL, ylab = NULL, par.settings = list(axis.line = list(col = 'black')),
                       panel = function(...) {
                           lattice::panel.levelplot(...)
                           lattice::panel.abline(h = start_coord, lty = 'dotted', col = 'black')
                           lattice::panel.abline(v = start_coord, lty = 'dotted', col = 'black')
                       })
}

#' Plot dendrogram
#'
#' @param htads `htad` object returned by `call_HTADs`.
#' @examples
#' load('data/chromosome18_10Mb.Rdata')
#' htads <- call_HTADs(chromosome18_10Mb)
#' plot_dendro(htads)
#' @export

plot_dendro <- function(htads) {
    plot(htads$dendro, labels = FALSE, hang = -1)
}

#' Plot retained variance
#'
#' @param input_data `data.frame` with 3 columns containing HiC data in the format `(bin1, bin2, score)`.
#' @param max_pcs The maximum number of principal components whose cumulative variance is to be plotted.
#' @param mark A place on the abcissa where to put a visual mark.
#' @examples
#' load('data/chromosome18_10Mb.Rdata')
#' plot_var(chromosome18_10Mb)
#' @export

plot_var <- function(input_data, max_pcs = NULL, mark = 200) {
  mat <- load_mat(input_data)
  
  # Sparse matrix and correlation.
  correlation_matrix <- sparse_cor(mat)$cor
  correlation_matrix[is.na(correlation_matrix)] <- 0 # Clean NA/NaN values.
  
  # PCA (compute first `number_pca` components).
  if (is.null(max_pcs)) max_pcs <- nrow(mat)
  number_pca <- min(max_pcs, nrow(mat))
  pca <- prcomp(correlation_matrix, rank. = number_pca)
  
  perc <- cumsum(pca$sdev^2) / sum(pca$sdev^2) * 100
  
  ggpubr::ggline(data.frame('variable' = 1:length(perc),
                            'value' = perc), 
         x = 'variable', 
         y = 'value', 
         plot_type = 'l',
         legend = 'none') +
    geom_segment(aes(x = mark, y = 0, xend = mark, yend = 100),
                 color = 'black', linetype = 'dotted', size = 0.1)
}

#' Call hierarchical TADs
#'
#' Computes a constrained hierarchical clustering of genomic regions in a HiC experiment,
#' choosing the optimal amount of information from the HiC matrix and selecting the most informative number of TADs.
#' @param input_data `data.frame` with 3 columns containing HiC data in the format `(bin1, bin2, score)`.
#' @param cores When `method` is `"accurate"`, the number of cores to use for parallel execution.
#' @param max_pcs The maximum number of principal components to retain for the analysis.
#' @param method Which version of the algorithm to use.
#' @param n_samples When `method` is `"fast"`, the number of samples used to approximate the optimal solution.
#' @param min_clusters ###### WRITE THIS UP ######.
#' @param plot Logical. Whether to plot the scores of every tested `n_pcs`/`n_clusters` combination.
#' @return `htad` object that defines the clustering of genomic regions.
#' @examples
#' load('data/chromosome18_10Mb.Rdata')
#' htads <- call_HTADs(chromosome18_10Mb)
#' @export

call_HTADs <- function(input_data, cores = 1, max_pcs = 200, method = c('accurate', 'fast'), n_samples = 60, min_clusters = 2) {
  # Load and clean data.
  mat <- load_mat(input_data)

  # Sparse matrix and correlation.
  # sparse_matrix <- Matrix::Matrix(mat , sparse = TRUE)
  # correlation_matrix <- sparse_cor(sparse_matrix)$cor
  correlation_matrix <- sparse_cor(mat)$cor
  correlation_matrix[is.na(correlation_matrix)] <- 0 # Clean NA/NaN values.

  # PCA (compute first `number_pca` components).
  number_pca <- min(max_pcs, nrow(mat))
  pca <- prcomp(correlation_matrix, rank. = number_pca)

  # Find optimal clustering parameters based on Calinhara score.
  method <- match.arg(method)
  optimal_params <- switch(method,
                           accurate = find_params_accurate(pca, number_pca, cores, min_clusters),
                           fast = find_params_fast(pca, number_pca, n_samples, min_clusters))

  # Cluster the PCs subset with the best mean-CHI criterion.
  pcs <- as.matrix(pca$x[, 1:optimal_params$n_PCs])
  row.names(pcs) <- 1:nrow(pcs)
  clust <- rioja::chclust(dist(pcs))

  # Create hierarchical cluster object.
  htads <- structure(list('n_pcs' = optimal_params$n_PCs[[1]],
                          'optimal_n_clusters' = optimal_params$n_clusters[[1]],
                          'dendro' = clust,
                          'clusters' = list(),
                          'scores' = optimal_params$scores),
                     class = 'htads')
  for (n in which(!is.na(optimal_params$scores[optimal_params$n_PCs, ]))) {
    eb <- cumsum(table(cutree(clust, k = n)))
    htads$clusters[[as.character(n)]] <- list('CH-index' = optimal_params$scores[optimal_params$n_PCs, n],
                                              'coord' = data.frame('start' = c(1, eb[-length(eb)] + 1, use.names = FALSE),
                                                                   'end' = eb))
  }

  htads
}
