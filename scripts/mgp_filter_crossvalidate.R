## Filter log-likelihoods from R phylopomp on simulated genealogies, for scripts/mgp_filter_crossvalidate.jl.
##
##   Rscript scripts/mgp_filter_crossvalidate.R <model> <ntrees> <Np> <nrep> <out.tsv> [lib]
##
## model is sir, si2r, lbdp, bdei or bdss. Each tree is simulated with runX, written with newick()
## and parsed back with parse_newick(), and the R pomp object is built from the parsed tree, so R and
## Julia read the same rounded string. Trees are kept when they have one root at t0 and 3 to 8 samples.
## Columns: model, tree, newick, time, nodetimes (sorted, incl. the end time), R logmeanexp, its SE,
## and lbdp_exact (NA for other models).
## The optional sixth argument is an R library that holds a different phylopomp build.
## SEED in the environment overrides the seed. Not part of test/runtests.jl: it needs R and phylopomp.
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]; ntrees <- as.integer(args[2]); Np <- as.integer(args[3])
nrep <- as.integer(args[4]); out <- args[5]
lib <- if (length(args) >= 6) args[6] else NULL
suppressMessages(library(phylopomp, lib.loc = lib))
suppressMessages(library(pomp))
set.seed(as.integer(Sys.getenv("SEED", "20261007")))

T <- c(sir = 4, si2r = 3, lbdp = 4, bdei = 4, bdss = 3)[[model]]
sim <- switch(model,
  sir  = function() runSIR(time = T, t0 = 0, Beta = 3, gamma = 1, psi = 0.3, omega = 0,
                           pop = 100, S0 = 0.99, I0 = 0.01, R0 = 0),
  si2r = function() runSI2R(time = T, t0 = 0, Beta = 4, kappa = 3, gamma = 1, omega = 0.5, chi = 0.3,
                            etaL = 0.5, etaH = 1, pop = 100, S0 = 0.99, IL0 = 0.01, IH0 = 0, R0 = 0),
  lbdp = function() runLBDP(time = T, t0 = 0, lambda = 1.5, mu = 0.5, psi = 0.3, chi = 0.2, n0 = 1),
  bdei = function() runBDEI(time = T, t0 = 0, sigma = 1, lambda = 2, mu = 0.5, chi = 0.4,
                            pop = 1, E0 = 0, I0 = 1),
  bdss = function() runBDSS(time = T, t0 = 0, lambda_nn = 1, lambda_ns = 0.3, lambda_sn = 1.5,
                            lambda_ss = 2.5, mu = 0.5, chi = 0.4, pop = 1, N0 = 1, S0 = 0),
  stop("unknown model: ", model))
mkpomp <- switch(model,
  sir  = function(y) sir_pomp(y, Beta = 3, gamma = 1, psi = 0.3, omega = 0,
                              S0 = 0.99, I0 = 0.01, R0 = 0, pop = 100),
  si2r = function(y) si2rs_pomp(y, Beta = 4, kappa = 3, gamma = 1, omega = 0.5, chi = 0.3,
                                etaL = 0.5, etaH = 1, S0 = 0.99, IL0 = 0.01, IH0 = 0, R0 = 0, pop = 100),
  lbdp = function(y) lbdp_pomp(y, lambda = 1.5, mu = 0.5, psi = 0.3, chi = 0.2, n0 = 1),
  bdei = function(y) bdei_pomp(y, sigma = 1, lambda = 2, mu = 0.5, chi = 0.4, pop = 1, E0 = 0, I0 = 1),
  bdss = function(y) bdss_pomp(y, lambda_nn = 1, lambda_ns = 0.3, lambda_sn = 1.5, lambda_ss = 2.5,
                               mu = 0.5, chi = 0.4, pop = 1, N0 = 1, S0 = 0))

## Simulate every tree first, so that the trees depend only on the seed and not on the filter.
trees <- list(); tries <- 0
while (length(trees) < ntrees && tries < 100000) {
  tries <- tries + 1
  s <- newick(sim())
  if (!nzchar(s)) next
  y <- parse_newick(s, t0 = 0, time = T)
  gi <- gendat(y)
  if (gi$nroot != 1 || gi$nsample < 3 || gi$nsample > 8 || gi$nodetime[1] != 0) next
  trees[[length(trees) + 1]] <- list(s = s, y = y, gi = gi)
}
rows <- character(0)
for (k in seq_along(trees)) {
  s <- trees[[k]]$s; y <- trees[[k]]$y; gi <- trees[[k]]$gi
  p <- mkpomp(y)
  ll <- replicate(nrep, logLik(pfilter(p, Np = Np)))
  est <- logmeanexp(ll, se = TRUE)
  ex <- if (model == "lbdp") lbdp_exact(y, lambda = 1.5, mu = 0.5, psi = 0.3, chi = 0.2, n0 = 1) else NA
  rows <- c(rows, paste(model, k, s, T, paste(sprintf("%.17g", sort(gi$nodetime)), collapse = ","),
                        sprintf("%.17g", est[1]), sprintf("%.17g", est[2]), sprintf("%.17g", ex), sep = "\t"))
}
k <- length(trees)
writeLines(rows, out)
cat(model, ": ", k, " trees from ", tries, " simulations\n", sep = "")
