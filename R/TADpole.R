# Imports:
# bigmemory
# Matrix
# doParallel
# dendextend (for plot)
# parallel
# foreach
# rioja
# fpc

#' Load a Hi-C matrix from a file
#'
#' @param mat_file path to the input file. Must be in a tab-delimited matrix format.
#' @param bad_frac fraction of the matrix to falg as bad rows/columns.
#' @param centromere_search split the matrix by the centrormere into two, smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) datasets.
#' @param hist_bad_columns plot the distribution of row/column coverage to help in selecting a useful value for `bad_frac`. Mostly for debugging.
#' @export

load_mat <- function(mat_file, bad_frac = 0.01, centromere_search = FALSE, hist_bad_columns = FALSE) {
    mat <- bigmemory::read.big.matrix(mat_file, type = 'double', sep = '\t')[, ]

    mat[is.na(mat)] <- 0 # Clean NA/NaN values.
    mat <- as.matrix(Matrix::forceSymmetric(mat, uplo = 'U'))
    rownames(mat) <- 1:nrow(mat)
    colnames(mat) <- 1:ncol(mat)

    # Detect bad columns.
    r <- rowMeans(mat)
    bad_columns <- diag(mat) == 0
    if (bad_frac) bad_columns <- bad_columns | r < quantile(r, seq(0, 1, by = bad_frac))[2]

    if (hist_bad_columns) hist(r, breaks = 50)

    attr(mat, 'bad_columns') <- names(which(bad_columns))
    message(paste(sum(bad_columns), 'bad columns found at position(s):'))
    message(paste(attr(mat, 'bad_columns'), collapse = ' '))

    if (centromere_search) {
        idx <- as.numeric(attr(mat, 'bad_columns'))

        list_bad_columns <- split(idx, cumsum(seq_along(idx) %in% (which(diff(idx) > 1) + 1)))
        centromere_start <- head(list_bad_columns[[which.max(lengths(list_bad_columns))]], 1)
        centromere_end <- tail(list_bad_columns[[which.max(lengths(list_bad_columns))]], 1)
        message(paste('centromere position:', centromere_start, centromere_end))
        if (centromere_start == 1 || centromere_end == nrow(mat)) {
            message('longest stretch of bad rows/columns at the ends, not splitting the matrix.')
            return(mat[!bad_columns, !bad_columns])
        }

        idx_p <- 1:(centromere_start - 1)
        idx_q <- (centromere_end + 1):nrow(mat)
        mat_p <- mat[idx_p, idx_p]
        mat_q <- mat[idx_q, idx_q]
        bad_colums_p <- idx[idx < centromere_start]
        bad_colums_q <- idx[idx > centromere_end]
        if (length(bad_colums_p)) mat_p <- mat_p[-bad_colums_p, -bad_colums_p]
        if (length(bad_colums_q)) mat_q <- mat_q[-bad_colums_q, -bad_colums_q]

        attr(mat_p, 'bad_columns') <- bad_colums_p
        attr(mat_q, 'bad_columns') <- bad_colums_q

        return(list(p = mat_p, q = mat_q, centromere = centromere_start:centromere_end))

    } else return(mat[!bad_columns, !bad_columns])
}

sparse_cor <- function(x) {
    # Create a sparse correlation matrix.
    covmat <- (crossprod(as.matrix(x)) - nrow(x) * tcrossprod(colMeans(x))) / (nrow(x) - 1)
    sdvec <- sqrt(diag(covmat))
    cormat <- covmat / tcrossprod(sdvec)
    list(cov = covmat, cor = cormat)
}

find_params <- function(pca, number_pca, min_clusters) {
    doParallel::registerDoParallel(cores = parallel::detectCores())
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

    # scores <- as.matrix(bigmemory::big.matrix(nrow = length(calinhara_score),
    #                                           ncol = max(sapply(calinhara_score, length)),
    #                                           type = 'integer',
    #                                           init = 0))
    scores <- matrix(nrow = length(calinhara_score), ncol = max(sapply(calinhara_score, length)))
    for (pc in 1:length(calinhara_score)) scores[pc, 1:length(calinhara_score[[pc]])] <- calinhara_score[[pc]]
    rownames(scores) <- 1:number_pca
    colnames(scores) <- 1:ncol(scores) # + 1

    optimal_PCs <- which.max(rowMeans(scores, na.rm = TRUE))
    optimal_n_clusters <- which.max(scores[optimal_PCs, ])
    message(paste('Optimal number of PCs:', optimal_PCs))
    message(paste('Optimal number of clusters:', optimal_n_clusters))

    list(n_PCs = optimal_PCs, n_clusters = optimal_n_clusters, scores = scores)
}

#' Plot dendrogram
#'
#' @param TADpole `TADpole` object returned as by function `TADpole`.
#' @examples
#' tadpole <- TADpole('data/chromosome18_10Mb.tsv')
#' plot_dendro(tadpole)
#' @export

plot_dendro <- function(tadpole) {
    dend <- as.dendrogram(tadpole$dendro)
    hpk <- dendextend::heights_per_k.dendrogram(dend)
    plot(cut(dend, h = hpk[tadpole$optimal_n_clusters])$upper, 
         leaflab = 'none',
         main="Dendrogram of all levels validated by the Broken-Stick model")
    # plot(cut(as.dendrogram(tadpole$dendro, hang = 10), h = tadpole$optimal_n_clusters)$upper, leaflab = 'none')
    rect.hclust(tadpole$dendro, k = tadpole$optimal_n_clusters)
    # cutted <- cut(dend, h = hpk[tadpole$optimal_n_clusters], ordered_result = FALSE)$upper
    # labels(cutted) <- 1:tadpole$optimal_n_clusters
}

#' Plot borders
#'
#' @param TADpole `TADpole` object returned as by function `TADpole`.
#' @param mat_file path to the input file. Must be in a tab-delimited matrix format.
#' @param centromere_search split the matrix by the centrormere into two, smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) datasets.
#' @examples
#' tadpole <- TADpole('data/chromosome18_10Mb.tsv')
#' plot_borders(tadpole, input_data = 'data/chromosome18_10Mb.tsv',centromere_search = FALSE)
#' @export

plot_borders <- function(tadpole, input_data, centromere_search) {
    mat <- read.big.matrix(input_data,type = "double",sep = "\t")
    mat <- as.matrix(input_data)
    mat[is.na(mat)] <- 0 # Clean NA/NaN values.
    mat = as.matrix(Matrix::forceSymmetric(mat, uplo = 'U'))
    colnames(mat) = seq(1:dim(mat)[1])
    row.names(mat) = seq(1:dim(mat)[1])
    print(paste("Dimension of the matrix:",dim(mat)))
    
    if (centromere_search == FALSE){
    
        start_coord <- tadpole$clusters[[as.character(tadpole$optimal_n_clusters)]]$coord$start
        end_coord <- tadpole$clusters[[as.character(tadpole$optimal_n_clusters)]]$coord$end}
    
    if (centromere_search == TRUE){
    
        start_coord <- tadpole$merging_arms$coord$start
        end_coord <- tadpole$merging_arms$coord$end}
    
    colors <- colorRampPalette(c('white', 'firebrick3'))
    
    lattice::levelplot(as.matrix(log(mat)), col.regions = colors, scales = list(draw = FALSE), colorkey = FALSE,
                       xlab = NULL, ylab = NULL, par.settings = list(axis.line = list(col = 'black')),
                       panel = function(...) {
                           lattice::panel.levelplot(...)
                           lattice::panel.abline(h = unique(c(start_coord - 0.5, end_coord + 0.5)), lty = 'dotted', col = 'black')
                           lattice::panel.abline(v = unique(c(start_coord - 0.5, end_coord + 0.5)), lty = 'dotted', col = 'black')
                       })}




#' Call hierarchical TADs
#'
#' Computes a constrained hierarchical clustering of genomic regions in a HiC experiment,
#' choosing the optimal amount of information from the HiC matrix and selecting the most informative number of TADs.
#' @param mat_file path to the input file. Must be in a tab-delimited matrix format.
#' @param max_pcs The maximum number of principal components to retain for the analysis.
#' @param min_clusters Minimum number of clusters into which partition the chromosome.
#' @param bad_frac fraction of the matrix to falg as bad rows/columns.
#' @param centromere_search split the matrix by the centrormere into two smaller matrices representing the chromosomal arms. Useful when woring with big (>15000 bins) datasets.
#' @param hist_bad_columns plot the distribution of row/column coverage to help in selecting a useful value for `bad_frac`. Mostly for debugging.
#' @return `tadpole` object that defines the clustering of genomic regions.
#' @details The `centromere_search` parameter will split the matrix into two by the region with the longes stretch of bad (low signal) rows/columns.
#' It will do so regardless of whether this stretch represents a true centromere or not. Note that this feature is useful when processing an entire chromosome,
#' but be cautious of interpreting the partitions as the two chromosomal arms (p and q) when working with smaller regions.
#' @examples
#' tadpole <- TADpole('data/chromosome18_10Mb.tsv')
#' @export

TADpole <- function(mat_file, max_pcs = 200, min_clusters = 2, bad_frac = 0.01, centromere_search = FALSE, hist_bad_columns = FALSE) {
    # Load and clean data.
    mat <- load_mat(mat_file, bad_frac = bad_frac, centromere_search = centromere_search, hist_bad_columns = hist_bad_columns)

    if (centromere_search) {
        fixed_clusters_arms <- c()
        names_clusters_arms <- c()
        tadpole <- structure(list(), class = 'tadpole')

        centromer <- mat$centromer
        for (arm in c('p', 'q')) {
            message(paste('Processing arm', arm))
            bad_columns <- attr(mat[[arm]], 'bad_columns')

            # Sparse matrix and correlation.
            correlation_matrix <- sparse_cor(mat[[arm]])$cor
            correlation_matrix[is.na(correlation_matrix)] <- 0 # Clean NA/NaN values.

            # PCA (compute first `number_pca` components).
            number_pca <- min(max_pcs, nrow(mat[[arm]]))
            pca <- prcomp(correlation_matrix, rank. = number_pca)

            # Find optimal clustering parameters based on Calinhara score.
            optimal_params <- find_params(pca, number_pca, min_clusters)

            # Cluster the PCs subset with the best mean-CHI criterion.
            pcs <- as.matrix(pca$x[, 1:optimal_params$n_PCs])
            clust <- rioja::chclust(dist(pcs))

            tadpole[[arm]]$n_pcs <- optimal_params$n_PCs
            tadpole[[arm]]$n_optimal <- optimal_params$n_clusters
            tadpole[[arm]]$dendro <- clust

            # Saving all the clusters accepted by the broken stick model
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

                tadpole[[arm]]$cluster[[as.character(k)]] <- list('coord' = coord)
            }

            # Merging optimal clusters in both arms.
            good_clusters <- cutree(clust, k = optimal_params$n_clusters)

            if (!is.null(bad_columns)) {
                # Add bad columns to the clusters.
                bad_clusters <- rep(0, length(bad_columns))
                names(bad_clusters) <- bad_columns

                # Merge bad columns with the original data.
                clusters <- c(good_clusters, bad_clusters)
                clusters <- clusters[order(as.numeric(names(clusters)))]

                rle_clusters <- fix_values(rle(clusters))
                fixed_clusters <- inverse.rle(rle_clusters)

            } else {
                # Merge bad columns with the original data.
                clusters <- good_clusters
                clusters <- clusters[order(as.numeric(names(clusters)))]

                rle_clusters <- fix_values(rle(clusters))
                fixed_clusters <- inverse.rle(rle_clusters)
            }

            fixed_clusters_arms <- c(fixed_clusters_arms, fixed_clusters, rep(0,length(centromer)))
            names_clusters_arms <- c(names_clusters_arms, names(clusters), as.character(centromer))
        }

        eb <- cumsum(rle(fixed_clusters_arms[1:(length(fixed_clusters_arms) - length(centromer))])$lengths)
        coord <- data.frame('start' = c(1, eb[-length(eb)] + 1, use.names = TRUE),
                            'end' = eb)
        coord <- coord[rle(fixed_clusters_arms[1:(length(fixed_clusters_arms) - length(centromer))])$values != 0, ]
        tadpole$merging_arms <- list('coord' = coord)

    } else {
        bad_columns <- attr(mat, 'bad_columns')

        # Sparse matrix and correlation.
        correlation_matrix <- sparse_cor(mat)$cor
        correlation_matrix[is.na(correlation_matrix)] <- 0 # Clean NA/NaN values.

        # PCA (compute first `number_pca` components).
        number_pca <- min(max_pcs, nrow(mat))
        pca <- prcomp(correlation_matrix, rank. = number_pca)

        # Find optimal clustering parameters based on Calinhara score.
        optimal_params <- find_params(pca, number_pca, min_clusters)

        # Cluster the PCs subset with the best mean-CHI criterion.
        pcs <- as.matrix(pca$x[, 1:optimal_params$n_PCs])
        clust <- rioja::chclust(dist(pcs))

        # Create hierarchical cluster object.
        tadpole <- structure(list('n_pcs' = optimal_params$n_PCs,
                                  'optimal_n_clusters' = optimal_params$n_clusters,
                                  'dendro' = clust,
                                  'clusters' = list(),
                                  'scores' = optimal_params$scores),
                             class = 'tadpole')

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

            tadpole$clusters[[as.character(k)]] <- list('CH-index' = optimal_params$scores[optimal_params$n_PCs, k],
                                                        'coord' = coord)
        }
    }

    tadpole
}

fix_values <- function(r) {
    zeros <- which(r$values == 0)
    zeros <- zeros[zeros != 1 & zeros != length(r$values)]
    for (i in zeros) {
        if (r$values[i - 1] == r$values[i + 1]) r$values[i] <- r$values[i - 1]
    }
    r
}
