## Milestone 3: distribution-level cross-validation of the Julia forward
## simulator (src/simulate.jl) against R phylopomp's runSEIR.
##
## Generates N independent SEIR genealogies from the R side and writes their
## Newick representations (one per line, empty line = total extinction, no
## samples) to a file for the Julia side (scripts/seir_crossvalidate.jl) to
## parse and compare against its own, independently-generated draws.
##
## Parameters are matched EXACTLY to scripts/seir_crossvalidate.jl's Julia
## side. In particular E0 = 0 (not R's default of 0.05) is deliberate: R's
## rinit (yaml/seir.yml) grafts round(pop*E0) lineages into Exposed and
## round(pop*I0) into Infectious; with E0 = 0 this reduces to a single
## founding lineage (into Infectious only), matching the single-root case
## the Julia simulator's acceptance tests already exercise. Using R's
## default E0 = 0.05 would require the (implemented but so far untested,
## see check_milestone1.md #6) multi-root/forest path on the Julia side --
## deliberately avoided here to keep this comparison isolated to the single
## engine question, not compounded with an unverified code path.
##
## Usage: Rscript scripts/seir_crossvalidate.R [N] [outfile]

suppressMessages(library(phylopomp))

args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args) >= 1) as.integer(args[1]) else 1000L
outfile <- if (length(args) >= 2) args[2] else "seir_crossvalidate_r_trees.txt"

Beta  <- 4.0
sigma <- 1.0
gamma <- 1.0
psi   <- 0.30
chi   <- 0.0
omega <- 1.0
pop   <- 100
S0    <- 0.99
E0    <- 0.0
I0    <- 0.01
R0    <- 0.0
t0    <- 0
time  <- 20.0

set.seed(20260813)

trees <- character(N)
for (i in seq_len(N)) {
  x <- runSEIR(
    time = time, t0 = t0,
    Beta = Beta, sigma = sigma, gamma = gamma, psi = psi, chi = chi, omega = omega,
    pop = pop, S0 = S0, E0 = E0, I0 = I0, R0 = R0
  )
  trees[i] <- newick(x)  ## "" (empty string) when the epidemic went extinct
}

writeLines(trees, outfile, sep = "\n")
cat(sprintf(
  "Wrote %d trees to %s (%d extinct/empty)\n",
  N, outfile, sum(nchar(trees) == 0)
))
