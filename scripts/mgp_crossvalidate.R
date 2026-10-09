## Reference draws from R phylopomp for scripts/mgp_crossvalidate.jl.
## One Newick string per line. A blank line means the run had no samples.
## <out>.states has the final population of each run.
## For MERS, SI2R, MTBD, BDEI and BDSS, the optional fourth file is
## `newick(x, obscure = FALSE)`, which keeps the deme of each sample.
## Parameters match the Julia script.
##
##   Rscript scripts/mgp_crossvalidate.R seir 1000 /tmp/r_seir.txt
##   Rscript scripts/mgp_crossvalidate.R seirchi 1000 /tmp/r_seirchi.txt
##   Rscript scripts/mgp_crossvalidate.R mers 2000 /tmp/r_mers.txt /tmp/r_mers_unobs.txt
##   Rscript scripts/mgp_crossvalidate.R sir 2000 /tmp/r_sir.txt
##   Rscript scripts/mgp_crossvalidate.R si2r 2000 /tmp/r_si2r.txt /tmp/r_si2r_unobs.txt
##   Rscript scripts/mgp_crossvalidate.R mtbd 2000 /tmp/r_mtbd.txt /tmp/r_mtbd_unobs.txt
##   Rscript scripts/mgp_crossvalidate.R lbdp 2000 /tmp/r_lbdp.txt
##   Rscript scripts/mgp_crossvalidate.R bdei 2000 /tmp/r_bdei.txt /tmp/r_bdei_unobs.txt
##   Rscript scripts/mgp_crossvalidate.R bdss 2000 /tmp/r_bdss.txt /tmp/r_bdss_unobs.txt
##
## Not part of test/runtests.jl. It needs R and phylopomp, and a KS test
## does not belong in CI.
suppressMessages(library(phylopomp))
args <- commandArgs(trailingOnly = TRUE)
model <- args[1]; N <- as.integer(args[2]); out <- args[3]
out2 <- if (length(args) >= 4) args[4] else NA
## SEED in the environment overrides the default seed
set.seed(as.integer(Sys.getenv("SEED", "20261003")))
trees <- character(N); trees2 <- character(N); states <- character(N)
## final population state from the yaml block "state:" ... "genealogy:",
## one "name=value ..." line per tree, written to <out>.states
final_state <- function(x) {
  ln <- strsplit(yaml(x), "\n")[[1]]
  a <- which(ln == "state:"); b <- which(ln == "genealogy:")
  kv <- trimws(ln[(a + 1):(b - 1)])
  paste(sub(": ", "=", kv), collapse = " ")
}
for (i in seq_len(N)) {
  if (model == "seir") {
    x <- runSEIR(time = 20, t0 = 0, Beta = 4, sigma = 1, gamma = 1, psi = 0.30, chi = 0,
                 omega = 1, pop = 100, S0 = 0.99, E0 = 0, I0 = 0.01, R0 = 0)
  } else if (model == "seirchi") {
    ## destructive sampling on: exercises the @mgp SEIR `culling` event
    x <- runSEIR(time = 20, t0 = 0, Beta = 4, sigma = 1, gamma = 1, psi = 0.20, chi = 0.10,
                 omega = 1, pop = 100, S0 = 0.99, E0 = 0, I0 = 0.01, R0 = 0)
  } else if (model == "mers") {
    ## rinit: Sc = round(Nc/(Sc0+Ic0)*Sc0) = 19, Ic = 1, Sh = 20, Ih = 0
    x <- runMERS(time = 10, t0 = 0,
                 Beta_cc = 3, Beta_ch = 0.5, Beta_hc = 0.5, Beta_hh = 3,
                 gamma_c = 1, gamma_h = 1, chi_c = 0.3, chi_h = 0.3,
                 Bc = 0.5, Bh = 0.5, Sc0 = 0.95, Sh0 = 1, Ic0 = 0.05, Ih0 = 0, Nc = 20, Nh = 20)
  } else if (model == "sir") {
    ## sampling keeps the host infectious
    x <- runSIR(time = 20, t0 = 0, Beta = 4, gamma = 1, psi = 0.3, omega = 0,
                pop = 100, S0 = 0.99, I0 = 0.01, R0 = 0)
  } else if (model == "si2r") {
    ## sampling removes the host (sample_death): PhyloPOMP.SI2R with r = 1
    x <- runSI2R(time = 20, t0 = 0, Beta = 4, kappa = 3, gamma = 1, omega = 0.5, chi = 0.3,
                 etaL = 0.5, etaH = 1, pop = 100, S0 = 0.99, IL0 = 0.01, IH0 = 0, R0 = 0)
  } else if (model == "mtbd") {
    ## r_1 < 1: some type-1 samples keep their lineage
    x <- runMTBD2(time = 6, t0 = 0, lambda_1_1 = 1.2, lambda_1_2 = 0.3, lambda_2_1 = 0.2,
                  lambda_2_2 = 0.9, m_1_2 = 0.2, m_2_1 = 0.1, mu_1 = 0.5, mu_2 = 0.5,
                  psi_1 = 0.3, psi_2 = 0.3, r_1 = 0.7, r_2 = 1, I1_0 = 1, I2_0 = 0)
  } else if (model == "lbdp") {
    x <- runLBDP(time = 4, t0 = 0, lambda = 1.5, mu = 0.5, psi = 0.3, chi = 0.2, n0 = 1)
  } else if (model == "bdei") {
    x <- runBDEI(time = 4, t0 = 0, sigma = 1, lambda = 2, mu = 0.5, chi = 0.4, pop = 1, E0 = 0, I0 = 1)
  } else if (model == "bdss") {
    x <- runBDSS(time = 3, t0 = 0, lambda_nn = 1, lambda_ns = 0.3, lambda_sn = 1.5, lambda_ss = 2.5,
                 mu = 0.5, chi = 0.4, pop = 1, N0 = 1, S0 = 0)
  } else stop("unknown model: ", model)
  trees[i] <- newick(x)                                      # prune=TRUE, obscure=TRUE
  if (!is.na(out2)) trees2[i] <- newick(x, obscure = FALSE)  # keeps sample demes
  states[i] <- final_state(x)
}
writeLines(trees, out)
writeLines(states, paste0(out, ".states"))
if (!is.na(out2)) writeLines(trees2, out2)
cat(sprintf("%s: wrote %d trees to %s (%d with no samples)\n", model, N, out, sum(nchar(trees) == 0)))
