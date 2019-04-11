estimate_resolution <- function(input_data) {
    resolution <- DescTools::GCD(input_data$V1, input_data$V2)
    print(paste('Estimated resolution:', resolution))
    resolution
}

#' Convert a 3-column data frame to a symmetric HiC matrix
#'
#' @param input_data `data.frame` with 3 columns containing HiC data in the format `(bin1, bin2, score)`.
#' @param empty_frac maximum fraction of empty bins allowed in a row/column.
#' @export

load_mat <- function(input_data, bad_frac = 0.01, plot_hist = FALSE, check = TRUE) {
    colnames(input_data) <- paste0('V', 1:3)

    # Build the full matrix by filling the missing coordinates.
    mat_coordinates <- unique(c(input_data$V1, input_data$V2))
    all_coordinates <- seq(min(mat_coordinates), max(mat_coordinates), estimate_resolution(input_data))
    full_df <- expand.grid(all_coordinates, all_coordinates)
    colnames(full_df) <- c('V1', 'V2')
    full_df$V3 <- 0
    mat <- xtabs(V3 ~ V1 + V2, aggregate(V3 ~ V1 + V2, rbind(input_data, full_df), sum))
    mat <- as.matrix(mat)
    mat[is.na(mat)] <- 0 # Clean NA/NaN values.

    # Check matrix symmetry.
    lower <- c(mat[lower.tri(mat)])
    upper <- c(mat[upper.tri(mat)])

    if (isSymmetric(unname(mat))) {
        print('Input matrix is already symmetric. Doing nothing.')
    } else if (!sum(upper)) {
        print('Filling the upper triangle of the matrix')
        mat = as.matrix(Matrix::forceSymmetric(mat, uplo = 'L'))
    } else if (!sum(lower)) {
        print('Filling the lower triangle of the matrix')
        mat = as.matrix(Matrix::forceSymmetric(mat, uplo = 'U'))
    } else if (check) stop('Input matrix is not symmetric!')

    # Detect bad columns.
    # bad_columns <- diag(mat) == 0 | rowMeans(mat == 0) > empty_frac)
    r <- rowMeans(mat)
    bad_columns <- diag(mat) == 0
    if (bad_frac) bad_columns <- bad_columns | r < quantile(r, seq(0, 1, by = bad_frac))[2]

    if (plot_hist) hist(r, breaks = 50)

    if (any(bad_columns)) {
        mat <- mat[!bad_columns, !bad_columns]
        attr(mat, 'bad_columns') <- names(which(bad_columns))

        print(paste(sum(bad_columns), 'bad columns found at position(s):'))
        print(attr(mat, 'bad_columns'))
    }

    mat
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
  calinhara_score <- foreach::`%dopar%`(foreach::foreach(i = 1:number_pca), {
    pcs <- as.matrix(pca$x[, 1:i])
    row.names(pcs) <- 1:nrow(pcs)

    clust <- rioja::chclust(dist(pcs))

    # Broken stick.
    bs <- rioja::bstick(clust, nrow(pcs) - 1, plot = FALSE)
    r <- rle(bs$dispersion > bs$bstick)
    n_cluster <- r$lengths[r$values][1]

    score <- rep(NA, n_cluster)
    min_clusters <- min(min_clusters, n_cluster)
    for (n in min_clusters:n_cluster) {
      chc <- cutree(clust, k = n)
      score[n] <- fpc::calinhara(pca$x, chc, cn = n)
    }

    score
  })

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

#' Save Bed file with the TAD's borders coordinates.
#'
#' @param htads `htad` object returned by `call_HTADs`.
#' @param chrom `character` with the name of the chromosome("chrX").
#' @param file path to save the file in bed format.
#' @examples
#' load('data/chromosome18_10Mb.Rdata')
#' htads <- call_HTADs(chromosome18_10Mb)
#' save_bed_coordinates(htads,chrom = "chr18","/scratch/User/chr18_borders.bed")
#' @export

save_bed_coordinates <- function(htads, chrom, output){
    start_coord <- htads$clusters[[as.character(htads$optimal_n_clusters)]]$coord$start
    end_coord <- htads$clusters[[as.character(htads$optimal_n_clusters)]]$coord$end
    bed_format <- cbind(chrom,start_coord,end_coord)
    write.table(bed_format,output,
                col.names = F,
                row.names = F,
                quote = F, sep = "\t")}


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
    mat <- load_mat(input_data, bad_frac = 0)
    start_coord <- htads$clusters[[as.character(htads$optimal_n_clusters)]]$coord$start
    end_coord <- htads$clusters[[as.character(htads$optimal_n_clusters)]]$coord$end
    colors <- colorRampPalette(c('white', 'firebrick3'))

    lattice::levelplot(as.matrix(log(mat)), col.regions = colors, scales = list(draw = FALSE), colorkey = FALSE,
                       xlab = NULL, ylab = NULL, par.settings = list(axis.line = list(col = 'black')),
                       panel = function(...) {
                           lattice::panel.levelplot(...)
                           lattice::panel.abline(h = unique(c(start_coord - 0.5, end_coord + 0.5)), lty = 'dotted', col = 'black')
                           lattice::panel.abline(v = unique(c(start_coord - 0.5, end_coord + 0.5)), lty = 'dotted', col = 'black')
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
    plot(cut(as.dendrogram(htads$dendro), h = htads$optimal_n_clusters)$upper, labels = FALSE, hang = -1)
    # rect.hclust(htads$dendro, k = htads$optimal_n_clusters)
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

plot_var <- function(input_data, max_pcs = NULL, mark = 200, percent = 0.8) {
  mat <- load_mat(input_data, percent)[[1]]

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
    ggplot2::geom_segment(ggplot2::aes(x = mark, y = 0, xend = mark, yend = 100),
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
#' @param min_clusters Minimum number of clusters into which partition the chromosome.
#' @param plot Logical. Whether to plot the scores of every tested `n_pcs`/`n_clusters` combination.
#' @return `htad` object that defines the clustering of genomic regions.
#' @examples
#' load('data/chromosome18_10Mb.Rdata')
#' htads <- call_HTADs(chromosome18_10Mb)
#' @export

call_HTADs <- function(input_data, cores = 1, max_pcs = 200, method = c('accurate', 'fast'), n_samples = 60, min_clusters = 1, bad_frac = 0.01, check = TRUE) {
  # Load and clean data.
  mat <- load_mat(input_data, bad_frac = bad_frac, check = check)
  bad_columns <- attr(mat, 'bad_columns')

  # Sparse matrix and correlation.
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
  clust <- rioja::chclust(dist(pcs))

  # Create hierarchical cluster object.
  htads <- structure(list('n_pcs' = optimal_params$n_PCs[[1]],
                          'optimal_n_clusters' = optimal_params$n_clusters[[1]],
                          'dendro' = clust,
                          'clusters' = list(),
                          'scores' = optimal_params$scores),
                     class = 'htads')

  for (k in which(!is.na(optimal_params$scores[optimal_params$n_PCs, ]))) {
      good_clusters <- cutree(clust, k = k)

      if (!is.null(bad_columns)) {
          # Add bad columns to the clusters.
          bad_clusters <- rep(0, length(bad_columns))
          names(bad_clusters) <- bad_columns

          # Merge bad columns with the original data.
          clusters <- c(good_clusters, bad_clusters)
          clusters <- clusters[order(as.numeric(names(clusters)))]

          rle_clusters <- fix_values(rle(clusters))
          fixed_clusters <- inverse.rle(rle_clusters)

          eb <- cumsum(rle(fixed_clusters)$length)
          coord <- data.frame('start' = c(1, eb[-length(eb)] + 1, use.names = FALSE),
                              'end' = eb)
          coord <- coord[rle(fixed_clusters)$values != 0, ]
      } else {
          eb <- cumsum(table(good_clusters))
          coord <- data.frame('start' = c(1, eb[-length(eb)] + 1, use.names = FALSE),
                              'end' = eb)
      }

      row.names(coord) <- NULL
      htads$clusters[[as.character(k)]] <- list('CH-index' = optimal_params$scores[optimal_params$n_PCs, k],
                                                'coord' = coord)
  }

  htads
}

fix_values <- function(r) {
    zeros <- which(r$values == 0)
    zeros <- zeros[zeros != 1 & zeros != length(r$values)]
    for (i in zeros) {
        if (r$values[i - 1] == r$values[i + 1]) r$values[i] <- r$values[i - 1]
    }
    r
}
