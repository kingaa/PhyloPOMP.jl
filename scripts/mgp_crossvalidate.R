## Cross-validation of the Julia forward simulator against R phylopomp:
## R-side reference draws for the SEIR and MERS models.
##
## Writes one Newick string per line (empty line = no samples). For MERS a
## second file with `obscure = FALSE` keeps the species (deme) of each sample
## node so the Julia side can compare camel/human sample counts.
##
## Parameters are matched EXACTLY to scripts/mgp_crossvalidate.jl.
##
## Usage:
##   Rscript scripts/mgp_crossvalidate.R seir 1000 /tmp/r_seir.txt
##   Rscript scripts/mgp_crossvalidate.R mers 2000 /tmp/r_mers.txt /tmp/r_mers_unobs.txt
##   julia --project=. scripts/mgp_crossvalidate.jl seir 1000 /tmp/r_seir.txt
##   julia --project=. scripts/mgp_crossvalidate.jl mers 2000 /tmp/r_mers.txt /tmp/r_mers_unobs.txt
##
## Deliberately not part of test/runtests.jl: needs R + phylopomp, and a
## stochastic KS test has no place in CI.
suppressMessages(library(phylopomp))
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]; N <- as.integer(args[2]); out <- args[3]
out2 <- if (length(args) >= 4) args[4] else NA
set.seed(20261003)
trees <- character(N); trees2 <- character(N)
for (i in seq_len(N)) {
  if (model == "seir") {
    ## chi = 0: the @mgp SEIR table has no destructive-sampling event.
    x <- runSEIR(time = 20, t0 = 0, Beta = 4, sigma = 1, gamma = 1, psi = 0.30, chi = 0,
                 omega = 1, pop = 100, S0 = 0.99, E0 = 0, I0 = 0.01, R0 = 0)
  } else if (model == "mers") {
    ## rinit: Sc = round(Nc/(Sc0+Ic0)*Sc0) = 19, Ic = 1, Sh = 20, Ih = 0
    x <- runMERS(time = 10, t0 = 0,
                 Beta_cc = 3, Beta_ch = 0.5, Beta_hc = 0.5, Beta_hh = 3,
                 gamma_c = 1, gamma_h = 1, chi_c = 0.3, chi_h = 0.3,
                 Bc = 0.5, Bh = 0.5, Sc0 = 0.95, Sh0 = 1, Ic0 = 0.05, Ih0 = 0, Nc = 20, Nh = 20)
  } else stop("unknown model: ", model)
  trees[i] <- newick(x)                                      # prune=TRUE, obscure=TRUE
  if (!is.na(out2)) trees2[i] <- newick(x, obscure = FALSE)  # keeps sample demes
}
writeLines(trees, out)
if (!is.na(out2)) writeLines(trees2, out2)
cat(sprintf("%s: wrote %d trees to %s (%d with no samples)\n", model, N, out, sum(nchar(trees) == 0)))
