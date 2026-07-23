#!/usr/bin/env Rscript
# Generate a handful of single-root SEIR trees for CBLV cross-check
# between phylodeep (Python) and PhyloPOMP.jl (Julia).
suppressPackageStartupMessages(library(phylopomp))

set.seed(77)
n_trees <- 5L
out_dir <- file.path("output", "crosscheck")
dir.create(file.path(out_dir, "trees"), recursive = TRUE, showWarnings = FALSE)

Beta  <- 4; sigma <- 1; gamma <- 1; psi <- 0
p_sample <- 0.4
chi   <- gamma * p_sample / (1 - p_sample)
pop   <- 1000; S0 <- 0.999; E0 <- 0; I0 <- 0.001; R0 <- 0

trees <- list()
for (i in seq_len(n_trees)) {
  ok <- FALSE
  attempts <- 0L
  while (!ok) {
    attempts <- attempts + 1L
    if (attempts > 200) stop("exceeded 200 attempts for tree ", i)
    t_end <- runif(1, 8, 25)
    x <- tryCatch(suppressWarnings(runSEIR(
      Beta = Beta, sigma = sigma, gamma = gamma,
      psi = psi, chi = chi, omega = 0,
      S0 = S0, E0 = E0, I0 = I0, R0 = R0,
      pop = pop, time = t_end
    )), error = function(e) NULL)
    if (is.null(x)) next
    info <- getInfo(x, newick = TRUE, nsample = TRUE, nroot = TRUE)
    if (info$nroot == 1 && info$nsample >= 20 && info$nsample <= 180) {
      ok <- TRUE
      nwk_file <- sprintf("trees/tree_%04d.nwk", i)
      writeLines(info$newick, file.path(out_dir, nwk_file))
      trees[[i]] <- data.frame(
        index = i, n_tips = info$nsample, n_root = info$nroot,
        nwk_file = nwk_file, p_sample = p_sample,
        stringsAsFactors = FALSE
      )
      cat(sprintf("tree %d: tips=%d attempts=%d\n", i, info$nsample, attempts))
    }
  }
}

manifest <- do.call(rbind, trees)
write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote %d trees to %s\n", n_trees, out_dir))
