# 1 bad col at 1 start
# 1 bad col at w middle
# 2 bad cols at both ends

library(data.table)

random_bed <- function(coords, size, bad_columns) {
    bins <- (1:size)[-bad_columns]
    borders <- sort(sample(bins[-1], nrow(coords) - 1))
    data.frame(coords[, 1],
               c(1, borders - 1),
               c(borders - 2, size))
}

tad_bin_index <- function(bed, size) {
    tad_wt <- rep(0, size)
    for (tad in 1:nrow(bed)) {
        for(bin in seq(bed[tad, 2], bed[tad, 3]) - bed[1, 2] + 1) {
            tad_wt[bin] <- tad
        }
    }
    tad_wt
}

scores <- function(tad_x, tad_y) {
    # TODO: Bad columns should not be counted.
    # stopifnot(length(tad_x) == length(tad_y))
    size <- length(tad_x)
    score <- 0
    scores <- c()
    bin1_all <- c()
    levels <- c()
    for(bin1 in 1:size) {
        x <- tad_x[bin1] != tad_x
        y <- tad_y[bin1] != tad_y
        score <- sum(xor(x, y))

        # print(c(bin1, score))
        scores = c(scores, score)
        bin1_all = c(bin1_all, bin1)
        levels = c(levels, i)

    }
    if (max(scores) == 0) cbind(bin1_all, scores, levels)
    else cbind(bin1_all, scores / max(scores), levels)
}

# Example.
start <- 293
end <- 486
size <- end - start + 1
bad_columns <- sort(sample(1:size, 3))
coords_wt <- read.table('../../../Downloads/9_level_wt.bed')
coords_inv1 <- read.table('../../../Downloads/9_level_inv1.bed')
index_wt <- tad_bin_index(coords_wt)
index_inv1 <- tad_bin_index(coords_inv1)
scores_diff_wt_inv1 <- scores(index_wt, index_inv1)

for (i in 1:1e4) {
        coord_random <- random_bed(coords_wt, size, bad_columns)
        index_random <- tad_bin_index(coord_random, size)
        difft_scores = scores(index_wt,index_random)
        write.table(difft_scores,
                    paste0("/scratch/pauli/TAD_caller/Capture_HiC_Kraft/ind1_random_overlapping/level_",i,"_",num,
                    ".bed"), col.names = F, row.names = F, quote = F)
}

write.table(scores_diff_wt_inv1,
            paste0("/scratch/pauli/TAD_caller/Capture_HiC_Kraft/ind1_overlapping/level_",i,".bed"),
                   col.names = F, row.names = F, quote = F)
