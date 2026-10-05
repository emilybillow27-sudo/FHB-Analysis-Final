# Regression checks for validation independence and identifier handling.
source("config.R")
source("scripts/functions.R")
set.seed(11)
markers <- matrix(sample(0:2, 60 * 25, replace = TRUE), 60, 25)
rownames(markers) <- paste0("G", seq_len(60))
grm <- sommer::A.mat(markers)
x <- data.frame(FullSampleName = rownames(grm),
                predicted.value = 100 + as.numeric(markers[, 1:3] %*% c(2, 3, 4)) + rnorm(60))
p <- predict_gblup(x[1:45, ], x[46:60, ], grm, "masking regression")
x[46:60, "predicted.value"] <- x[46:60, "predicted.value"] + 1000
q <- predict_gblup(x[1:45, ], x[46:60, ], grm, "changed withheld observations")
stopifnot(nrow(p) == 15L, all(is.finite(p$GBLUP)), mean(p$GBLUP) > 90,
          isTRUE(all.equal(p$GBLUP, q$GBLUP)))
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), "error"))
expect_error(predict_gblup(x[1:45, ], x[45:60, ], grm, "overlap"))
expect_error(prediction_input(data.frame(FullSampleName = c("G1", "G1"), TRAIT = "INC",
                                        predicted.value = c(1, 2)), "INC", grm))
ambiguous <- diag(2)
rownames(ambiguous) <- colnames(ambiguous) <- c("G-1", "G1")
expect_error(harmonize(data.frame(FullSampleName = "G1"), ambiguous))
stopifnot(is.na(safe_cor(c(1, 1, 1), c(1, 2, 3))))
cat("Validation independence, prediction scale, overlap, and ID checks passed.\n")
