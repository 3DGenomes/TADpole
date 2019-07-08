bin_index <- function(bed, size) {
    tad_index <- rep(0, size)
    for (tad in 1:nrow(bed)) {
        for(bin in seq(bed[tad, 2], bed[tad, 3]) - bed[1, 2] + 1) {
            tad_index[bin] <- tad
        }
    }
    tad_index
}

random_bed <- function(coords, size, bad_columns) {
    bins <- (1:size)[-bad_columns]
    borders <- sort(sample(bins[-1], nrow(coords) - 1))
    data.frame(coords[, 1],
               c(1, borders - 1),
               c(borders - 2, size))
}

#' Compute diffT score between two TAD calls
#'
#' @param `bed_x, bed_y` the two calls to compare. Each must be a `data.frame` with a BED-like format.
#' @export

diffT <- function(bed_x, bed_y) {
    if (nrow(bed_x) != nrow(bed_y)) stop('Both calls must have the same number of TADs.')

    start_x <- bed_x[1, 2]
    start_y <- bed_y[1, 2]
    end_x <- bed_x[nrow(bed_x), 3]
    end_y <- bed_y[nrow(bed_y), 3]

    tad_x <- bin_index(bed_x, end_x - start_x + 1)
    tad_y <- bin_index(bed_y, end_y - start_y + 1)

    # Extend terminal TADs for the missing bad columns at the beginning/end of the indices.
    tad_x <- c(rep(1, max(0, start_x - start_y)),
               tad_x,
               rep(max(tad_x), max(0, end_y - end_x)))
    tad_y <- c(rep(1, max(0, start_y - start_x)),
               tad_y,
               rep(max(tad_y), max(0, end_x - end_y)))
    # As a result, the indices of both calls should have the same length.
    stopifnot(length(tad_x) == length(tad_y))
    
    # Bad columns are counted as not matching.
    scores <- c()
    for(bin in 1:length(tad_x)) {
        x <- tad_x[bin] != tad_x | tad_x[bin] == 0
        y <- tad_y[bin] != tad_y | tad_y[bin] == 0
        scores <- c(scores, sum(xor(x, y)))
    }
    score_sum <- cumsum(scores)
    if (max(scores) == 0) score_sum
    else score_sum / max(score_sum)
}
