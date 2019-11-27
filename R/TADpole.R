#' Load a Hi-C matrix from a file
#'
#' @param mat_file path to the input file. Must be in a tab-delimited matrix format.
#' @param chr string with the chromosome name.
#' @param start numeric start position of the region or of the chromosome.
#' @param end numeric end position of the region or of the chromosome.
#' @param resol numeric resolution/binning of the Hi-C experiment.
#' @param bad_frac fraction of the matrix to falg as bad rows/columns.
#' @param centromere_search split the matrix by the centrormere into two, smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) datasets.
#' @examples
#' mat_file <- system.file("extdata", "raw_chr18:460-606_20kb.tsv", package = "TADpole")
#' mat <- load_mat(mat_file, chr = "chr18", start = 496, end = 606, resol = 20000)
#' @export

load_mat <- function(mat_file, chr, start, end, resol, bad_frac = 0.01, centromere_search = FALSE) {

    mat <- as.matrix(bigmemory::read.big.matrix(mat_file, type = 'double', sep = '\t')[, ])

    mat[is.na(mat)] <- 0 # Clean NA/NaN values.
    mat <- as.matrix(Matrix::forceSymmetric(mat, uplo = 'U'))
    rownames(mat) <- 1:nrow(mat)
    colnames(mat) <- 1:ncol(mat)

    # Plot Hi-C map
    colors <- colorRampPalette(c("white", "firebrick3"))

    p1 = (lattice::levelplot(as.matrix(log(mat)),
                   main=list(paste0("Raw Hi-C contact map \n",chr,":",start,"-",end),
                             side=1,line=0.2,fontsize=11),
                   sub= list(paste0("Resolution:",resol),fontsize=9),
                   col.regions = colors, scales = list(draw = FALSE), colorkey = TRUE,
                   xlab = NULL, ylab = NULL, par.settings = list(axis.line = list(col = 'black'))))

    # Detect bad columns.
    r <- rowMeans(mat)
    bad_columns <- diag(mat) == 0
    if (bad_frac) bad_columns <- bad_columns | r < quantile(r, seq(0, 1, by = bad_frac))[2]

    r_melt = reshape2::melt(r)
    gghist <- ggpubr::gghistogram(r_melt,  x = "value", rug = TRUE,bins = 50,add_density = TRUE,
                  fill = "#00AFBB") +
    ggplot2::geom_vline(xintercept = (quantile(r, seq(0, 1, by = bad_frac))[[2]]),
           linetype = 3)

     p2 = (ggpubr::ggpar(gghist,submain = "... Bad columns",
     main = "Interactions counts",
     font.main = c(11, "bold"),
     font.submain = c(9, "bold"),
     font.x = c(8, "bold"),
     font.y = c(8,"bold"),
     xlab = "Frequency of Hi-C interactions", ylab = "Counts"))

    gridExtra::grid.arrange(p1, p2, ncol=2, heights=c(3.5,2), widths=c(2,1))

    message(paste(sum(bad_columns), 'bad columns found at position(s):'))
    message(paste(names(which(bad_columns)), collapse = ' '))

    if (any(bad_columns) && centromere_search == TRUE) {
        attr(mat, 'bad_columns') <- names(which(bad_columns))
        idx <- as.numeric(attr(mat, 'bad_columns'))

        list_bad_columns <- split(idx, cumsum(seq_along(idx) %in% (which(diff(idx) > 1) + 1)))
        centromere_start <- head(list_bad_columns[[which.max(lengths(list_bad_columns))]], 1)
        centromere_end <- tail(list_bad_columns[[which.max(lengths(list_bad_columns))]], 1)
        message(paste('centromere position:', centromere_start, centromere_end))
        if (centromere_start == 1 || centromere_end == nrow(mat)) {
            message('longest stretch of bad rows/columns at the ends, not splitting the matrix.')
            mat <- mat[!bad_columns, !bad_columns]
            attr(mat, 'bad_columns') <- names(which(bad_columns))
            return(mat)
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

    } else {
        mat <- mat[!bad_columns, !bad_columns]
        attr(mat, 'bad_columns') <- names(which(bad_columns))
        return(mat)
    }
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

    scores <- matrix(nrow = length(calinhara_score), ncol = max(sapply(calinhara_score, length)))
     # scores <- as.matrix(bigmemmory::big.matrix(nrow = length(calinhara_score),
     #                             ncol = max(sapply(calinhara_score, length)),
     #                             type = "integer", init = 0))

    for (pc in 1:length(calinhara_score)) scores[pc, 1:length(calinhara_score[[pc]])] <- calinhara_score[[pc]]
    rownames(scores) <- 1:number_pca
    colnames(scores) <- 1:ncol(scores) # + 1

    optimal_PCs <- which.max(rowMeans(scores, na.rm = TRUE))
    optimal_n_clusters <- which.max(scores[optimal_PCs, ])
    message(paste('Optimal number of PCs:', optimal_PCs))
    message(paste('Optimal number of clusters:', optimal_n_clusters))

    list(n_PCs = optimal_PCs, n_clusters = optimal_n_clusters, scores = scores)
}

#' Plot hierarchical_plot
#' @param mat_file path to the input file. Must be in a tab-delimited matrix format.
#' @param TADpole `TADpole` object returned as by function `TADpole`.
#' @param chr string with the chromosome name.
#' @param start numeric start position of the region or of the chromosome.
#' @param end numeric end position of the region or of the chromosome.
#' @param resol numeric resolution/binning of the Hi-C experiment.
#' @param centromere_search split the matrix by the centrormere into two, smaller matrices representing the chromosomal arms. Useful when working with big (>15000 bins) datasets.

#' @examples
#' mat_file <- system.file("extdata", "raw_chr18:460-606_20kb.tsv", package = "TADpole")
#' tadpole <- TADpole(mat_file, chr = "chr18", start = 496, end = 606, resol = 20000)
#' plot_hierarchy(mat_file, tadpole, chr = "chr18", start = 496, end = 606, resol = 20000)
#' @export

plot_hierarchy <- function(mat_file, tadpole, chr, start, end, resol, centromere_search = FALSE) {

    # Matrix partition
    mat <- as.matrix(bigmemory::read.big.matrix(mat_file, type = 'double', sep = '\t')[, ])
    mat[is.na(mat)] <- 0 # Clean NA/NaN values.
    mat <- as.matrix(Matrix::forceSymmetric(mat, uplo = 'U'))
    rownames(mat) <- 1:nrow(mat)
    colnames(mat) <- 1:ncol(mat)

    colors <- colorRampPalette(c('white', 'firebrick3'))

    if (centromere_search) {
        start_coord <- tadpole$merging_arms$start
        end_coord <- tadpole$merging_arms$end
    } else {
        start_coord <- tadpole$clusters[[as.character(tadpole$optimal_n_clusters)]]$start
        end_coord <- tadpole$clusters[[as.character(tadpole$optimal_n_clusters)]]$end
    }

    matrix_partition <- lattice::levelplot(as.matrix(log(mat)),
                  # main=list('Hierarchical chromatin organization\nchr18:9,200,000-12,130,000',side=1,line=0.5),
                  # sub= paste0("Optimal number of PCs:",as.character(tadpole$n_pcs),"       ",
                   #            "Optimal number of clusters:",as.character(tadpole$optimal_n_clusters),"\n",
                    #           "____ Optimal partition             .... Significant partitions"),
                   col.regions = colors, scales = list(draw = FALSE), colorkey = TRUE,
                   xlab = NULL, ylab = NULL, par.settings = list(axis.line = list(col = 'black')),
                   panel = function(...) {
                       lattice::panel.levelplot(...)
                       for (i in seq(length(tadpole$clusters))){
                       lattice::panel.points(tadpole$clusters[[as.character(i)]]$start - 0.5,
                                             tadpole$clusters[[as.character(i)]]$end + 0.5,
                                             col="black",
                                             type="s", cex=4, lty="dashed")}
                       lattice::panel.points(start_coord - 0.5,
                                             end_coord + 0.5,
                                             col="blue",
                                             type="s", cex=4)})

    if (centromere_search) {
      ggdraw(xlim = c(0, 1.3), ylim = c(0, 1.3)) +
      draw_plot(matrix_partition, x = 0.2, y = 0.1, width = 0.9, height = 0.9) +
      draw_label(paste0('Hierarchical chromatin organization\n',chr,":",start,"-",end),
                 colour = "#80404080", size = 17,x = 0.7, y = 1.15) +
      draw_label(paste0("Optimal number of PCs in Q arm:",as.character(tadpole$q$n_pcs),"    ",
                        "Optimal number of PCs in P arm:",as.character(tadpole$p$n_pcs),"     ",'\n',
                        "Optimal number of clusters in Q arm:",as.character(tadpole$q$optimal_n_clusters), "   ",
                        "Optimal number of clusters in P arm:",as.character(tadpole$p$optimal_n_clusters),"\n",
                         "____ Optimal partition             .... Significant partitions"),
                 colour = "#80404080", size = 10,x = 0.7, y = 1)
    } else {
        cut <- length(tadpole$clusters)    # Number of clusters
        dendr <- ggdendro::dendro_data(tadpole$dendro, type="rectangle")
        clust <- cutree(tadpole$dendro, k = cut)               # find 'cut' clusters
        clust.df <- data.frame(label = names(clust), cluster = clust)

        # Split dendrogram into upper grey section and lower coloured section
        height <- unique(dendr$segments$y)[order(unique(dendr$segments$y), decreasing = TRUE)]
        cut.height <- mean(c(height[cut], height[cut-1]))
        dendr$segments$line <- ifelse(dendr$segments$y == dendr$segments$yend &
           dendr$segments$y > cut.height, 1, 2)
        dendr$segments$line <- ifelse(dendr$segments$yend  > cut.height, 1, dendr$segments$line)

        # Number the clusters
        dendr$segments$cluster <- c(-1, diff(dendr$segments$line))
        change <- which(dendr$segments$cluster == 1)
        for (i in 1:cut) dendr$segments$cluster[change[i]] = i + 1
        dendr$segments$cluster <-  ifelse(dendr$segments$line == 1, 1,
                     ifelse(dendr$segments$cluster == 0, NA, dendr$segments$cluster))
        dendr$segments$cluster <- zoo::na.locf(dendr$segments$cluster)

        # Consistent numbering between segment$cluster and label$cluster
        clust.df$label <- factor(clust.df$label, levels = levels(dendr$labels$label))
        clust.df <- plyr::arrange(clust.df, label)
        clust.df$cluster <- factor((clust.df$cluster), levels = unique(clust.df$cluster), labels = (1:cut) + 1)
        dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")

        # Positions for cluster labels
        n.rle <- rle(dendr$segments$cluster)
        N <- cumsum(n.rle$lengths)
        N <- N[seq(1, length(N), 2)] + 1
        N.df <- dendr$segments[N, ]
        N.df$cluster <- N.df$cluster - 1

        # Plot the dendrogram
        p <- ggplot2::ggplot() +
          ggplot2::geom_segment(data = ggdendro::segment(dendr),
                                ggplot2::aes(x=x, y=y, xend=xend, yend=yend, size=factor(line), colour=factor(cluster)),
              lineend = "square", show.legend = FALSE) +
           #scale_colour_manual(values = c("black", rainbow(cut))) +
          ggplot2::scale_size_manual(values = c(.1, 0.3)) +
          ggplot2::geom_text(data = N.df, ggplot2::aes(x = x, y = y, label = factor(cluster),  colour = factor(cluster + 1)),
              hjust = 1.5, show.legend = FALSE) +
          ggplot2::scale_y_reverse(expand = c(0.2, 0)) +
          ggplot2::labs(x = NULL, y = NULL) +
          ggplot2::coord_flip() +
          ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                axis.text.y = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank(),
                panel.background = ggplot2::element_rect(fill = "white"),
                panel.grid = ggplot2::element_blank())

        cowplot::ggdraw(xlim = c(0, 1.3), ylim = c(0, 1.3)) +
          cowplot::draw_plot(matrix_partition, x = 0.4, y = 0.1, width = 0.9, height = 0.9) +
          cowplot::draw_plot(p, x = 0, y = 0.08, width = .4, height = 0.9) +
          cowplot::draw_label(paste0('Hierarchical chromatin organization\n',chr,":",start,"-",end),
                     colour = "#80404080", size = 17,x = 0.7, y = 1.15) +
          cowplot::draw_label(paste0("Optimal number of PCs:",as.character(tadpole$n_pcs),"       ",
                            "Optimal number of clusters:",as.character(tadpole$optimal_n_clusters),"\n",
                             "____ Optimal partition             .... Significant partitions"),
                     colour = "#80404080", size = 11,x = 0.7, y = 1)
    }
}

#' Plot matrix of Calinski-Harabasz (CH) index
#' @param TADpole `TADpole` object returned as by function `TADpole`.

#' @examples
#' mat_file <- system.file("extdata", "raw_chr18:460-606_20kb.tsv", package = "TADpole")
#' tadpole <- TADpole(mat_file)
#' CH_map(tadpole)
#' @export

CH_map <- function(tadpole) {
  # TODO: make it compatible with centromere = TRUE.
    df = data.frame(Var2 = tadpole$n_pcs, Var1 = tadpole$optimal_n_clusters)
    s <- t(tadpole$scores)
    tadpole_melt <- reshape2::melt(s)
    tadpole_melt <-tadpole_melt[tadpole_melt$value!=0,]
    print(head(tadpole_melt))

    ggplot2::ggplot(tadpole_melt, ggplot2::aes(x = Var2, y = Var1)) +
      ggplot2::geom_raster(ggplot2::aes(fill=value)) +
      viridis::scale_fill_viridis() +
      ggplot2::labs(x="Number of PCs", y="Number of clusters", title='Caliski-Harabasz index') +
      ggplot2::theme_bw() + ggplot2::theme(axis.text.x=ggplot2::element_text(size=9, angle=0, vjust=0.3),
                           axis.text.y=ggplot2::element_text(size=9),
                           plot.title=ggplot2::element_text(size=11)) +
      ggplot2::geom_point(data = df,color = "blue",size = 1.5) +
      ggplot2::geom_vline(xintercept=tadpole$n_pcs, linetype="dashed", color = "blue",size = 0.5)
}

#' Call hierarchical TADs
#'
#' Computes a constrained hierarchical clustering of genomic regions in a HiC experiment,
#' choosing the optimal amount of information from the HiC matrix and selecting the most informative number of TADs.
#' @param mat_file path to the input file. Must be in a tab-delimited matrix format.
#' @param max_pcs The maximum number of principal components to retain for the analysis.
#' @param min_clusters Minimum number of clusters into which partition the chromosome.
#' @param bad_frac fraction of the matrix to falg as bad rows/columns.
#' @param chr string with the chromosome name.
#' @param start numeric start position of the region or of the chromosome.
#' @param end numeric end position of the region or of the chromosome.
#' @param resol numeric resolution/binning of the Hi-C experiment.
#' @param centromere_search split the matrix by the centrormere into two smaller matrices representing the chromosomal arms. Useful when woring with big (>15000 bins) datasets.
#' @return `tadpole` object that defines the clustering of genomic regions.
#' @details The `centromere_search` parameter will split the matrix into two by the region with the longes stretch of bad (low signal) rows/columns.
#' It will do so regardless of whether this stretch represents a true centromere or not. Note that this feature is useful when processing an entire chromosome,
#' but be cautious of interpreting the partitions as the two chromosomal arms (p and q) when working with smaller regions.
#' @examples
#' mat_file <- system.file("inst/extdata", "raw_chr18:460-606_20kb.tsv", package = "TADpole")
#' tadpole <- TADpole(mat_file, chr = "chr18", start = 496, end = 606, resol = 20000)
#' @export

TADpole <- function(mat_file, max_pcs = 200, min_clusters = 2, bad_frac = 0.01,
                    chr, start, end, resol, centromere_search = FALSE) {

    # Load and clean data.
    mat <- load_mat(mat_file, chr, start, end, resol,
                    bad_frac = bad_frac, centromere_search = centromere_search)

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
            tadpole[[arm]]$optimal_n_clusters <- optimal_params$n_clusters
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

                tadpole[[arm]]$cluster[[as.character(k)]] <- coord
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
        tadpole$merging_arms <- coord

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

            tadpole$clusters[[as.character(k)]] <- coord
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
