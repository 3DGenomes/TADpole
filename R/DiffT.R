bin_index <- function(bed, size) {
    tad_index <- rep(0, size)
    for (tad in 1:nrow(bed)) {
        for(bin in seq(bed[tad, 2], bed[tad, 3]) - bed[1, 2] + 1) {
            tad_index[bin] <- tad
        }
    }
    tad_index
}

#' Compute diffT score between two TAD calls
#' @param `bed_x,bed_y` two `data.frame`s with a BED-like format with 3 columns: chromosome, start and end coordinates of each TAD, in bins.
#' @examples
#' difft_control_case <- diffT(control, case)
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

#' Compute a random set of TAD border coordinates from a sample partition
#'
#' @param `bed` a `data.frame` with a BED-like format with 3 columns: chromosome, start and end coordinates of each TAD, in bins.
#' @param `bad_columns` a numeric `vector` with the positions of the bad columns. TAD borders will not be placed on bad columns. Default value of `NULL` means no bad columns will be introduced.
#' @examples
#' random_coords <- random_bed(control)
#' @export

random_bed <- function(bed, bad_columns = NULL) {
    start <- bed[1, 2]
    end <- bed[nrow(bed), 3]
    size <- end - start + 1

    if (is.null(bad_columns)) bins <- start:end
    else bins <- (start:end)[-bad_columns]

    borders <- sort(sample(bins[-1], nrow(bed) - 1))
    data.frame(chrom = bed[, 1],
               start = c(start, borders - 1),
               end = c(borders - 2, start + size - 1))
}
