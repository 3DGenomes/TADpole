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

    # Add 0s for the missing bad columns at the beginning/end of the indices.
    tad_x <- c(rep(0, max(0, start_x - start_y)),
               tad_x,
               rep(0, max(0, end_y - end_x)))
    tad_y <- c(rep(0, max(0, start_y - start_x)),
               tad_y,
               rep(0, max(0, end_x - end_y)))
    # As a result, the indices of both calls should have the same length.
    stopifnot(length(tad_x) == length(tad_y))

    # Bad columns are counted as not matching.
    size <- length(tad_x)
    scores <- c()
    for(bin in 1:size) {
        x <- tad_x[bin] != tad_x | tad_x[bin] == 0
        y <- tad_y[bin] != tad_y | tad_y[bin] == 0
        scores <- c(scores, sum(xor(x, y)))
    }
    if (max(scores) == 0) cbind(1:size, scores)
    else cbind(1:size, scores / max(scores))
}

# start <- 293
# end <- 486
# size <- end - start + 1
bad_columns <- sort(sample(1:size, 3))
coords_wt <- read.table('../../../Downloads/9_level_wt.bed')
coords_inv1 <- read.table('../../../Downloads/9_level_inv1.bed')

index_wt <- tad_bin_index(coords_wt, size)
index_inv1 <- tad_bin_index(coords_inv1, size)
diff_wt_inv1 <- diffT(index_wt, index_inv1)

coord_random <- random_bed(coords_wt, size, bad_columns)
index_random <- tad_bin_index(coord_random, size)
difft_scores <- diffT(index_wt,index_random)

for (num in seq(1,10,1))
        {
        htads_random_coord = (random_bed(i,htads_wt_coord))
        index_array_random = tad_bin_index(htads_random_coord,size)
        scores_diff = scores(index_array_wt,index_array_random)
     #   write.table(scores_diff,
      #              paste0("/scratch/pauli/TAD_caller/Capture_HiC_Kraft/ind4_random_overlapping/level_",i,"_",num,
      #              ".bed"), col.names = F, row.names = F, quote = F)
        }

i = 5

htads_wt_coord = fread(paste0("/scratch/pauli/TAD_caller/Capture_HiC_Kraft/wt_levels/",i,"_level.bed"))
htads_inv1_coord = fread(paste0("/scratch/pauli/TAD_caller/Capture_HiC_Kraft/inv1_levels/",i,"_level.bed"))

bins_clusters = inverse.rle(rle(c(htads_wt_coord$V3-htads_wt_coord$V2) + 1))
bins = rep(1:length(bins_clusters),bins_clusters)
print(length(bins))

